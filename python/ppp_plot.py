#!/usr/bin/env python3

import sys
import argparse


t_bias = 10  # I'd prefer to reduce the timers initially to 0, but sometimes slaves start before master


class PPP_Node:
    """A single node on tree of events with all event enties and links"""

    parent = None
    proc = None
    label = ""

    def __init__(self, label, time=None):
        self.label = label
        self.start = []
        self.stop = []
        if time is not None:
            self._set_time(time)
        self.children = {}  # this flattens the tree a bit but we can recontruct it from timings anyway, if needed
        self.parent = None

    def _set_time(self, time):
        if time > 0:
            if (len(self.start) != len(self.stop)):
                sys.stderr.write("Warning: orphaned start for '" + self.path() + "' " + str(time) + " vs " + str(self.start[-1]) + "\n")
                self.stop.append(None)
            self.start.append(time)
        elif time < 0:
            if (len(self.start) != len(self.stop) + 1):
                sys.stderr.write("Warning: orphaned stop for '" + self.path() + "' " + str(time) + " vs " + str(self.stop[-1]) + "\n")
                self.start.append(None)
            self.stop.append(-time)
        else:
            sys.stderr.write("Warning: time == 0 for '" + self.path() + "'\n")

    def _add(self, label, time):
        if time > 0:  # create/update child
            if label not in self.children:
                self.children[label] = PPP_Node(label, time)
                self.children[label].parent = self
            else:
                self.children[label]._set_time(time)
            return self.children[label]
        else:
            if label == self.label:
                self._set_time(time)
            else:
                sys.stderr.write("Warning: unfinished for '" + self.path() + "' " + str(time) + "\n")
            return self.parent  # most likely won't recover properly

    def path(self):
        return str(self.proc) if self.parent is None else self.parent.path() + "/" + self.label

    def shortpath(self):
        return "" if self.parent is None else self.parent.shortpath() + "/" + self.label

    def print(self, indent=0):  # need to filter through parent range
        try:
            for _ in range(min(len(self.start), len(self.stop))):
                print("  " * indent + "'" + self.label + "' %.6f %.6f" % (self.start[_] - t_bias, self.stop[_] - self.start[_]))
        except TypeError:
            print("  " * indent + "'" + self.label + "' TypeError: ", self.start, self.stop)
        for i in self.children:
            self.children[i].print(indent=indent + 1)

    def get_all_ev(self):
        evl = [[self.shortpath(), self.start, self.stop, self.own_time, self.child_time]] if len(self.start) > 0 else []
        for i in self.children:
            evl += self.children[i].get_all_ev()
        return evl

    def calculate_time(self):
        self.own_time = 0.
        for i in range(len(self.start)):
            self.own_time += self.stop[i] - self.start[i]
        self.child_time = 0.
        for _ in self.children:
            self.child_time += self.children[_].calculate_time()
        return self.own_time

    def apply_exclusions(self):
        if args.exclude is not None:
            for _ in list(self.children.keys()):
                if self.children[_].label in args.exclude:
                    del self.children[_]

        if args.maxdepth is not None:
            for _ in list(self.children.keys()):
                if self.shortpath().count("/") >= args.maxdepth:
                    del self.children[_]

        for _ in self.children:
            self.children[_].apply_exclusions()

    def root_matches(self):
        if self.label in args.root:
            return [self]
        else:
            match = []
            for _ in self.children:
                match += self.children[_].root_matches()
            return match

    def shift_time(self, dt):
        for i in range(len(self.start)):
            self.start[i] -= dt
            self.stop[i] -= dt
        for _ in self.children:
            self.children[_].shift_time(dt)


class PPP_Tree:
    """An event tree for one Piernik process"""

    def __init__(self, name):
        self.thread = name
        self.root = PPP_Node("/")  # root node is special and does not have time
        self.root.proc = self.thread
        self._last = self.root

    def _add(self, label, time):
        self._last = self._last._add(label, time)

    def print(self):
        print("Process: " + str(self.thread))
        self.root.print()

    def rebase_root(self):
        self.newroot = PPP_Node("/")
        self.newroot.proc = self.thread
        for r in self.root.root_matches():
            nn = r.parent.path() if r != self.root else "/"
            newname = nn
            i = 1
            while newname in self.newroot.children:  # find some unique label
                newname = nn + " %d" % i
                i += 1
            self.newroot.children[newname] = r
        self.root = self.newroot


class PPP:
    """A collection of event trees from one or many Piernik processes from a single run"""

    def __init__(self, name):
        self.name = name
        self.trees = {}  # separate event tree for each process

    def print(self):
        print(self.name)
        for i in sorted(self.trees):
            self.trees[i].print()

    def _decode_text(self):
        self.nthr = 0
        self.bigbang = None
        file = open(self.name, 'r')
        for line in file:
            l = line.split()
            if self.bigbang is None:
                self.bigbang = float(l[1]) - t_bias
            self._add(int(l[0]), line[line.index(l[2]):].strip(), float(l[1]) + self.bigbang * (-1. if float(l[1]) > 0. else 1.))  # proc, label, time
            if int(l[0]) >= self.nthr:
                self.nthr = int(l[0]) + 1
        self.descr = "'%s' (%d thread%s)" % (self.name, self.nthr, "s" if self.nthr > 1 else "")

    def _add(self, proc, label, time):
        if proc not in self.trees:
            self.trees[proc] = PPP_Tree(proc)
        self.trees[proc]._add(label, time)

    def calculate_time(self):
        self.total_time = 0.
        for _ in self.trees:
            self.trees[_].root.calculate_time()
            self.total_time += self.trees[_].root.child_time  # root does not have time

    def apply_exclusions(self):
        if args.exclude is not None or args.maxdepth is not None:
            for p in self.trees:
                self.trees[p].root.apply_exclusions()

    def rebase_root(self):
        if args.root is not None:
            earliest = None
            for p in self.trees:
                self.trees[p].rebase_root()
                for i in self.trees[p].root.children:
                    if earliest is None:
                        earliest = self.trees[p].root.children[i].start[0]
                    else:
                        earliest = min(earliest, self.trees[p].root.children[i].start[0])
            if earliest is not None:
                for p in self.trees:
                    self.trees[p].root.shift_time(earliest - t_bias)


class PPPset:
    """A collection of event trees from one or many Piernik runs"""

    def __init__(self, fnamelist):
        """Let's focus on ASCII decoding for a while. Don't bother with HDF5 until we implement such type of PPP dump"""
        self.run = {}
        for fname in fnamelist:
            self.run[fname] = PPP(fname)
        for _ in self.run:
            self.run[_]._decode_text()
            self.run[_].apply_exclusions()
            self.run[_].rebase_root()
            self.run[_].calculate_time()

    def print(self, otype):
        self.out = ""
        if otype == "tree":
            for _ in self.run:
                print("")
                self.run[_].print()
        elif otype == "gnu":
            self.print_gnuplot()
            if args.output is None:
                print(self.out)
            else:
                ofile = open(args.output[0], "w")
                ofile.write(self.out)
                ofile.close()
        elif otype == "summary":
            # ToDo: add own time, %time and "called from"
            print("ARGS: ", args)
            for f in self.run:
                ed = {}
                print("\n## File: '%s', %d threads, bigbang = %.7f" % (self.run[f].name, self.run[f].nthr, self.run[f].bigbang))
                for p in self.run[f].trees:
                    for e in self.run[f].trees[p].root.get_all_ev():
                        e_base = e[0].split('/')[-1]
                        if e_base not in ed:
                            ed[e_base] = [0, 0., 0.]  # calls, cumul. own time, cumul. child. time
                        ed[e_base][0] += len(e[1])
                        ed[e_base][1] += e[3]
                        ed[e_base][2] += e[4]
                ml = len(max(ed, key=len)) if len(ed) > 0 else 0
                print("# %-*s %20s %8s %15s \033[97m%20s\033[0m %10s" % (ml - 2, "label", "avg. CPU time (s)", "%time", "avg. occurred", "avg. non-child", "% of total"))
                print("# %-*s %20s %8s %15s \033[97m%20s\033[0m %10s" % (ml - 2, "", "(per thread)", "in child", "(per thread)", "time (s)", "time"))
                skipped = [0, 0.]
                for e in sorted(ed.items(), key=lambda x: x[1][1] - x[1][2], reverse=True):
                    if (e[1][1] - e[1][2]) / self.run[f].total_time > args.cutsmall[0] / 100.:
                        print("%-*s %20.6f %8s %15d%s \033[97m%20.6f\033[0m %10.2f" % (ml, e[0], e[1][1] / self.run[f].nthr,
                                                                                       ("" if e[1][2] == 0. else "%8.2f" % ((100 * e[1][2] / e[1][1]) if e[1][1] > 0. else 0.)),
                                                                                       e[1][0] / self.run[f].nthr,
                                                                                       " " if e[1][0] % self.run[f].nthr == 0 else "+",
                                                                                       (e[1][1] - e[1][2]) / self.run[f].nthr, 100. * (e[1][1] - e[1][2]) / self.run[f].total_time))
                    else:
                        skipped[0] += 1
                        skipped[1] += (e[1][1] - e[1][2]) / self.run[f].nthr
                if skipped[0] > 0:
                    print("# (skipped %d timers that contributed %.6f seconds of non-child time = %.2f%% of total time)" % (skipped[0], skipped[1], 100. * self.run[f].nthr * skipped[1] / self.run[f].total_time))

    def print_gnuplot(self):
        self.descr = ""
        ev = {}
        peff = 0
        gcnt = 0
        gomit = 0
        for f in self.run:
            for p in self.run[f].trees:
                evlist = self.run[f].trees[p].root.get_all_ev()
                for e in evlist:
                    e_base = e[0].split('/')[-1].split()[0]
                    if e_base not in ev:
                        ev[e_base] = {}
                    if p + peff not in ev[e_base]:
                        ev[e_base][p + peff] = []
                    if e[3] > args.cutsmall[0] / 100. * self.run[f].total_time / self.run[f].nthr:
                        ev[e_base][p + peff].append(e)
            peff = peff + self.run[f].nthr + (1 if len(self.run) > 1 else 0)  # spacing only when we have multiple files
            self.descr = self.run[f].descr + ("\\n" if len(self.descr) > 0 else "") + self.descr

        ndel = []
        for e in ev:
            n = 0
            for p in ev[e]:
                n += len(ev[e][p])
            if n == 0:
                ndel.append(e)
        for i in ndel:
            del(ev[i])
        if len(ndel) > 0:
            self.descr += "\\n(omitted %d timer%s below %g%% threshold)" % (len(ndel), "" if len(ndel) == 1 else "s", args.cutsmall[0])

        self.out += "#!/usr/bin/gnuplot\n$PPPdata << EOD\n\n"
        i_s = 1
        i_e = 2
        ev_cnt = {}
        for e in ev:
            ev_cnt[e] = 0
            for p in ev[e]:
                for t in range(len(ev[e][p])):
                    ev_cnt[e] += len(ev[e][p][t][i_s])
        for e in sorted(ev_cnt, key=lambda x: ev_cnt[x]):  # if we need to skip something due to MAXOUPUT, let it be the most abundant thing
            ftname = None
            for p in ev[e]:
                if ftname is None and len(ev[e][p]) > 0:
                    ftname = ev[e][p][0][0]
            self.out += "# __" + e + "__ '" + ftname + "'\n"
            for p in ev[e]:
                for t in range(len(ev[e][p])):
                    depth = ev[e][p][t][0].count("/") - 1
                    for _ in range(min(len(ev[e][p][t][i_s]), len(ev[e][p][t][i_e]))):
                        try:
                            if (ev[e][p][t][i_e][_] - ev[e][p][t][i_s][_]) > 0.:
                                if gcnt < args.maxoutput[0]:
                                    self.out += "%.7f %.7f %.7f %.7f %.7f %.7f %d %d\n" % (ev[e][p][t][i_s][_] - t_bias, depth,
                                                                                           ev[e][p][t][i_s][_] - t_bias, ev[e][p][t][i_e][_] - t_bias,
                                                                                           depth + float(p) / peff, depth + float(p + 1) / peff, depth, p)
                                else:
                                    gomit += 1
                                gcnt += 1
                        except TypeError:
                            sys.stderr.write("Warning: inclomplete event '", e, "' @", p, " #", str(t), ev[e][p][t])
                    self.out += "\n"
            self.out += "\n"
        if gomit > 0:
            self.descr += "\\n(omitted %d entries above maxoutput limit of %d)" % (gomit, args.maxoutput[0])
        self.out += "EOD\n\n# Suggested gnuplot commands:\nset key outside horizontal\n"
        self.out += "set xlabel 'time (walltime seconds)'\nset ylabel 'timer depth + proc/nproc'\n"
        self.out += 'set title "%s"\n' % self.descr.replace('_', "\\\\_")
        if len(ev) > 0:
            pline = "plot $PPPdata "
            nocomma = True
            fs = "solid 0.1"
            i = 0
            gp_colors = 8
            for e in ev:
                ind = e.split("/")[-1]
                if nocomma:
                    nocomma = False
                else:
                    pline += ', "" '
                pline += ' index "__' + ind + '__" with boxxy title "' + ind.replace('_', "\\\\\\_") + '" fs  ' + fs
                i += 1
                if i >= gp_colors:  # gnuplot seems to have only gp_colors colors for solid filling
                    npat = int(i / gp_colors)
                    if npat >= 3:  # pattern #3 does not show borders
                        npat += 1
                    fs = "pat %d" % npat
            self.out += pline + "\n"
            if len(ndel) > 0:
                self.out += 'print "Timers below threshold: ' + " ".join(ndel) + '"\n'
            self.out += "pause mouse close\n"
        else:
            self.out += 'show title; print "No data to plot"'

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description="Piernik Precise Profiling Presenter", epilog="""
examples:

plot profile of 64th step from file.ascii with gnuplot (hopefully in the interactive mode):
    ppp_plot.py file.ascii -r 'step 64'| gnuplot
    ppp_plot.py file.ascii -r 'step 64'-o file.gnu; gnuplot file.gnu

when gnuplot fails to set up desired teminal by default, try to set $GNUTERM (qt or x11 are recommended):
    ppp_plot.py file.ascii | GNUTERM=qt gnuplot

print list of top-lefel timers (steps) present in file.ascii:
    ppp_plot.py file.ascii -t -d 1
(the atep names are followed by their time offset and length)

find most time-consuming timers:
    ppp_plot.py file.ascii -s
(this also helps to identify most often called timers to be excluded in case of performance problems with gnuplot)

plot profile without the identified too-often called timer (e.g. "Loechner_mark") that makes interactive gnuplot to choke:
    ppp_plot.py file.ascii -e Loechner_mark | gnuplot
same as above but don't filter out timers that are contributing less than 0.1%%:
    ppp_plot.py file.ascii -e Loechner_mark -% 0 | gnuplot
""")

parser.add_argument("filename", nargs='+', help="PPP ascii file(s) to process")
parser.add_argument("-o", "--output", nargs=1, help="processed output file (gnuplot only)")
parser.add_argument("-%", "--cutsmall", nargs=1, default=[.1], type=float, help="skip contributions below CUTSMALL%% (default = 0.1%%)")
parser.add_argument("-e", "--exclude", nargs='+', help="do not show EXCLUDEd timer(s)")  # multiple excudes
parser.add_argument("-r", "--root", nargs='+', help="show only ROOT and their children")
parser.add_argument("-d", "--maxdepth", type=int, help="limit output to MAXDEPTH")
parser.add_argument("-m", "--maxoutput", nargs=1, default=[50000], type=int, help="limit output to MAXOUTPUT enries (gnuplot only, default = 50000)")
# parser.add_argument("-c", "--check", help="do a formal check only")

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-g", "--gnuplot", action="store_const", dest="otype", const="gnu", default="gnu", help="gnuplot output (default)")
pgroup.add_argument("-t", "--tree", action="store_const", dest="otype", const="tree", help="tree output")
pgroup.add_argument("-s", "--summary", action="store_const", dest="otype", const="summary", help="short summary")

args = parser.parse_args()
all_events = PPPset(args.filename)
all_events.print(args.otype)

# for j  in *.ppprofile.ascii ; do for i in `awk '{print $2}' $j |  sort |  uniq` ; do in=`grep "$i *  -" $j | wc -l` ; out=`grep "$i *  [1-9]" $j | wc -l` ; [ $in != $out ] && echo $j $i $in $out ; done |  column -t ; done
