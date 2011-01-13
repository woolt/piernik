#!/usr/bin/perl
# $Id$
#
# Usage: makemake {<program name> {<F90 compiler or fc or f77 or cc or c>}}
#
# Generate a Makefile from the sources in the current directory.  The source
# files may be in either C, FORTRAN 77, Fortran 90 or some combination of
# these languages.  If the F90 compiler specified is cray or parasoft, then
# the Makefile generated will conform to the conventions of these compilers.
# To run makemake, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# Written by Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
# Modified by Kacper Kowalik <kowalik@astri.umk.pl> October, 24, 2007
# Modified by Dominik Woltanski <minikwolt@astri.umk.pl> April-June, 2008
# Modified by Oskar Karczewski <olk@ncac.torun.pl> December, 2008

use Switch;
use File::Path;
use File::Copy;
use File::Basename;
use Getopt::Long;

my $setup_line = "./setup";

my %moo=();
my $OK = 0;
my $FAILED = 1;
my @ARGV_SAVED = @ARGV;

open DEF_CONF, "< ./.setuprc";
do {
  $_ = <DEF_CONF>;
  s/\n//;
  if ($_) { @ARGV = split(/ /); }
  else    { @ARGV = @ARGV_SAVED; }
  $setup_line .= " @ARGV";
  GetOptions("help"=>\$moo{help},
             "problems"=>\$moo{problems},
             "units"=>\$moo{units},
#             "dumpenv"=>\$moo{dumpenv},
             "define=s@"=>\@define,
             "param:s"=>\$moo{param},
             "compiler:s"=>\$moo{compiler},
             "obj:s"=>\$moo{obj},
             "linkexe"=>\$moo{linkexe},
             "verbose"=>\$moo{verbose},
             "copy"=>\$moo{copy},
             "nocompile"=>\$moo{nocompile},
             "laconic"=>\$moo{laconic});
} while ($_);
close DEF_CONF;

print "Called '$setup_line'\n\n";

if ($moo{verbose}) {
  print "Setup options:\n";
  foreach (keys %moo) { print "$_ = $moo{$_}\n"; }
  print "\n";
}

my $SELF = fileparse("$0");
my $argc = @ARGV;

# DEBUG
#print "Unprocessed by Getopt::Long\n" if $ARGV[0];
#foreach (@ARGV) {
#  print "$_\n";
#}

if ($argc < 1) {
  $moo{help} = 1;
  $is_prob = 0;
}

if ($moo{problems}) {
   system ( "cat ./problems/*/info");
   print " \n";
   exit $OK;
} elsif ($moo{units}) {
   system ( "grep uses ./src/base/constants.F90");
   print " \n";
   exit $OK;
} elsif ($moo{help}) {
   print "USAGE: " . $SELF . " [PROBLEM] [OPTIONS]\n";
   print "Available options:\n";
   print " --help            - print this help message and exit\n";
   print " --problems        - print available problems and exit\n";
   print " --units           - print available unit systems and exit\n";
   print " --copy            - hard-copy source files instead of linking them\n";
   print " --linkexe         - do not copy obj/piernik to runs/<problem> but link it instead\n";
#   print " --dumpenv\t- do not compile the source but dump the environment to \'obj/version.F90\'\n";
   print " --verbose         - try to confuse the user with some diagnostics ;-)\n";
   print " --nocompile       - Create object directory, check for circular dependencies, but do not compile or prepare run directory.
                     In combination with --copy will prepare portable object directory.\n";
   print " --laconic         - Compress stdout from make process.\n";
   print " --param <file>    - use <file> instead problem.par\n";
   print " --define <flag>   - inject precompiler directive to piernik.def\n";
   print " --obj <postfix>   - use obj_<postfix> directory instead of obj/ and runs/<problem>_<postfix>/
                     rather than runs/<problem>/\n";
   print " --compiler <*.in> - choose specified config from compilers directory\n\n";
   print "The file .setuprc may contain some preferred options to be included everytime,\nfor example '--linkexe --compiler <local_configuration> --laconic'\n";
   print "\n";

   print "all-in-one wrapper for a selected problem:\n";
   print "> " . $SELF . " <problem>\n";
   print "this will:\n";
   print "* create directory \'./obj\' with source files for <problem>\n";
   print "* compile the source using the default compiler\n";
   print "* put the executable in \'./runs/<problem>\'\n\n";

   print "to change the compiler settings (after \'" . $SELF . " <problem>\'):\n";
   print "> cd obj\n";
   print "> ./newcompiler <settingsname>\n";
   print "> make\n";
   print "then copy files to your run directory (optional), e.g.\n";
   print "> cp {piernik,problem.par} ../runs/<problem>\n";
   print "> cd ../runs/<problem>\n\n";

   print "to run PIERNIK:\n";
   print "edit problem.par as appropriate, e.g.\n";
   print "* add var names for visualisation => var(<number>)='<name>'\n";
   print "* change domain dimensions/resolution => DOMAIN_SIZES\n";
   print "* change domain divisions for parallel processing => MPI_BLOCKS\n";
   print "* change frequency of data dumps => dt_* entries\n";
   print "* etc.\n";
   print "execute\n";
   print "> ./piernik\n";
   print "or for <np> parallel processes\n";
   print "> mpirun -n <np> ./piernik\n\n";

   print "HEALTH WARNINGS:\n";
   print "* the contents of \'./obj\' and \'./runs/<problem>\' are overwritten\n";
   print "  each time \'" . $SELF . " <problem>\' is run, unless you specify -obj <postfix>\n";
   print "  in which case the contents of runs/<problem>_<postfix> will be only updated\n";
   print "* the def file \'piernik.def\' is copied only for reference, to change flags\n";
   print "  with which the source is compiled edit \'./problems/<problem>/piernik.def\'\n";
   print "* by default PIERNIK will read the configuration file \'problem.par\' from the\n";
   print "  working directory, to use alternative configurations execute\n";
   print "  \'./piernik <directory with an alternative problem.par>\'\n";
#   print "!===================================================================================!\n";
   print "\nEnjoy your work with the Piernik Code!\n";

   if ( $is_prob == 0 ) { exit $OK; }
   else                 { exit $FAILED; }
}

# defaults
if ($moo{compiler}) { $compiler = $moo{compiler}.".in"; }
else                { $compiler = 'default.in'; }

if ($moo{param}) { $param_file = $moo{param}; }
else             { $param_file = 'problem.par'; }

if ($moo{obj}) { $objdir = "obj_$moo{obj}"; }
else           { $objdir = "obj"; }

open MAKEIN, "< compilers/".$compiler or die $!;
@makein = <MAKEIN>;
close MAKEIN;

rmtree( $objdir );
mkpath( $objdir );
$probdir = "problems/" . $ARGV[0] . "/";
if(!-r $probdir) {
   print "$probdir: $!\n";
   print "problem \'" . $ARGV[0] . "\' does not exist\n";
   exit $FAILED;
}

@prob = <$probdir*.F90 $probdir$param_file>;
for (@prob) {
   s/problems/\.\.\/problems/;
}
@base = <src/base/*.F90 src/base/*.c src/base/*.h>;
for (@base) {
   s/src/\.\.\/src/;
}
@addons = ();
copy("compilers/newcompiler", "$objdir/newcompiler");
system("chmod a+x $objdir/newcompiler");
$defs = $probdir."piernik.def";
$odefs = "$objdir/piernik.def";
if ($moo{copy}) { copy($defs, $odefs) or die "failed copying piernik.def"; }
else {
  symlink("../".$defs, $odefs) or die "failed linking piernik.def";
  unless ( -e $odefs ) { die "$odefs: $!"; }
}
if(@define) {
   foreach $inc (@define) {
      push(@d, "".$inc."\n");
      $cppflags = $cppflags . " -D" . $inc;
   }
}
@d = grep (/define/, `echo \'#include \"$defs\"\' > foo.f90 && cpp $cppflags -dM foo.f90 && rm foo.*`);
if ($moo{verbose}) { print "Defined symbols:\n", @d, "\n"; }

if( grep { /GRAV/ }  @d) {
   push(@addons, "../src/gravity/gravity.F90");
   push(@addons, "../src/gravity/hydrostatic.F90");
}
if( grep { /POISSON_FFT/ }  @d) {push(@addons, "../src/gravity/poissonsolver.F90");}
if( grep { /PGPLOT/ }  @d) {push(@addons, "../src/vizualization/viz_pgplot.F90");}
if( grep { /MULTIGRID/ } @d) {
   @temp=<src/multigrid/multigrid*.F90>;
   if( grep { /COSM_RAYS/} @d) {push(@temp, <src/fluids/cosmicrays/multigrid*.F90>)};
   if( grep { /GRAV/ }  @d) {push(@temp, <src/gravity/multigrid*.F90>)};
   for (@temp) { s/src/\.\.\/src/;}
   push(@addons, @temp);
}
if( grep { /SHEAR/ }  @d) {push(@addons, "../src/shear/shear.F90");}
push(@addons, "../src/hdf5/dataio_hdf5.F90");
push(@addons, "../src/hdf5/list_hdf5.F90");
push(@addons, "../src/fluids/initfluids.F90");
push(@addons, "../src/fluids/fluidindex.F90");
push(@addons, "../src/fluids/fluxes.F90");
push(@addons, "../src/fluids/common/timestepfuncs.F90");

if( grep { /MAGNETIC/} @d) {
   push(@addons, "../src/magnetic/magboundaries.F90");
#   if ( grep { /UMUSCL/} @d) {} else {
      push(@addons, "../src/scheme/rtvd_split/advects.F90");
#   }
   if( grep { /RESIST/} @d) {push(@addons, "../src/magnetic/resist/resistivity.F90");} }

if( grep { /IONIZED/} @d) {
   push(@addons, "../src/fluids/ionized/initionized.F90");
   push(@addons, "../src/fluids/ionized/fluxionized.F90");
   push(@addons, "../src/fluids/ionized/timestepionized.F90"); }
if( grep { /NEUTRAL/} @d) {
   push(@addons, "../src/fluids/neutral/initneutral.F90");
   push(@addons, "../src/fluids/neutral/fluxneutral.F90");
   push(@addons, "../src/fluids/neutral/timestepneutral.F90"); }
if( grep { /DUST/} @d) {
   push(@addons, "../src/fluids/dust/initdust.F90");
   push(@addons, "../src/fluids/dust/fluxdust.F90");
   push(@addons, "../src/fluids/dust/timestepdust.F90"); }
if( grep { /COSM_RAYS/} @d) {
   push(@addons, "../src/fluids/cosmicrays/initcosmicrays.F90");
   push(@addons, "../src/fluids/cosmicrays/fluxcosmicrays.F90");
   push(@addons, "../src/fluids/cosmicrays/timestepcosmicrays.F90");
   push(@addons, "../src/fluids/cosmicrays/crdiffusion.F90");
   push(@addons, "../src/fluids/cosmicrays/crhelpers.F90");
   }
if( grep { /COSM_RAYS_SOURCES/} @d) {
   push(@addons, "../src/fluids/cosmicrays/sourcecosmicrays.F90");
   push(@addons, "../src/fluids/cosmicrays/crcomposition.F90");
   }
if( grep { /FLUID_INTERACTIONS/} @d) {
   push(@addons, "../src/fluids/interactions/interactions.F90");
   push(@addons, "../src/fluids/interactions/timestepinteractions.F90");
}
if( grep { /SN_SRC/} @d) {
   push(@addons, "../src/supernovae/snsources.F90");
}

#if( grep { /UMUSCL/} @d) {
#   push(@addons, "../src/scheme/unsplit_muscl_hancock/fluidupdate.F90");
#} else {
   push(@addons, "../src/scheme/rtvd_split/fluidupdate.F90");
   push(@addons, "../src/scheme/rtvd_split/sweeps.F90");
   push(@addons, "../src/scheme/rtvd_split/rtvd.F90");
#}

@files = ( @base, @prob, @addons );
@symln = ();
foreach $file (@files) {
   $pos = rindex($file,"\/")+1;
   $len = length($file);
   push(@symln, "$objdir/".substr($file,$pos,$len-$pos) );
}

for $i (0 .. $#symln){
  if ($moo{copy}) { copy   ("$objdir/".$files[$i], $objdir) or die "failed copying $files[$i]: $!"; }
  else {
    symlink($files[$i], $symln[$i]) or die "failed linking $files[$i]: $!";
    unless ( -e $symln[$i] ) { die "$symln[$i]: $!"; }
  }
}
chdir $objdir;
open(MAKEFILE, "> Makefile-prep");

# Source listing
#
#print MAKEFILE @makein;
print MAKEFILE "RM ?= /bin/rm\n";
print MAKEFILE "MV ?= /bin/mv\n";
print MAKEFILE "ifdef CHECK_MAGIC\n";
print MAKEFILE "\tRM  = /bin/true\n";
print MAKEFILE "\tMV  = /bin/true\n";
print MAKEFILE "\tF90 = /bin/true\n";
print MAKEFILE "endif\n";
print MAKEFILE "SRCS =\t";
@srcs = <*.c *.F90 version.F90>;
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Object listing
#
print MAKEFILE "OBJS =\t";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.o/ };
&PrintWords(8, 0, @objs);
print MAKEFILE "\n\n";
#
# Define common macros
#
print MAKEFILE "LIBS +=\${STATIC} -lhdf5_fortran -lhdf5hl_fortran -lz \${DYNAMIC}\n";
if( grep { /PGPLOT/ } @d)  {print MAKEFILE "LIBS += -lpgplot\n";}
if( grep { /SHEAR/ || /MULTIGRID/ } @d ) {print MAKEFILE "LIBS += `pkg-config --libs fftw3`\n";}
if( grep { /POISSON_FFT/ }  @d) {print MAKEFILE "LIBS += `pkg-config --libs fftw3` `pkg-config --libs lapack`\n";}
print MAKEFILE "CPPFLAGS += " . $cppflags . "\n";
#
# make
#
if ($moo{laconic}) { print MAKEFILE "SILENT = 1\n\n"; }
else               { print MAKEFILE "SILENT = 0\n\n"; }
print MAKEFILE "ifndef PRECOMP\n";
print MAKEFILE "\tPRECOMP=cpp\n";
print MAKEFILE "endif\n";
print MAKEFILE "ifeq (\"\$(SILENT)\",\"1\")\n";
print MAKEFILE "MAKEFLAGS += -s\n";
print MAKEFILE "define ECHO_FC\n";
print MAKEFILE "\@echo [FC] \$<\n";
print MAKEFILE "endef\n";
print MAKEFILE "define ECHO_CC\n";
print MAKEFILE "\@echo [CC] \$<\n";
print MAKEFILE "endef\n";
print MAKEFILE "else\n";
print MAKEFILE "define ECHO_FC\n";
print MAKEFILE "endef\n";
print MAKEFILE "define ECHO_CC\n";
print MAKEFILE "endef\n";
print MAKEFILE "endif\n\n";
print MAKEFILE ".PHONY: print_fc\n\n";
print MAKEFILE "all: date print_fc \$(PROG) \n\n";
print MAKEFILE "print_fc:\n";
print MAKEFILE "ifeq (\"\$(SILENT)\",\"1\")\n";
print MAKEFILE "\t\@echo FC = \$(F90) \$(CPPFLAGS) \$(F90FLAGS) -c\n";
print MAKEFILE "\t\@echo CC = \$(CC) \$(CPPFLAGS) \$(CFLAGS) -c\n";
print MAKEFILE "endif\n\n";
print MAKEFILE "\$(PROG): \$(OBJS)\n";
print MAKEFILE "\t\@echo \$(", &LanguageCompiler($ARGV[1], @srcs),") \$(LDFLAGS) -o \$@ '*.o' \$(LIBS)\n";
print MAKEFILE "\t\@\$(", &LanguageCompiler($ARGV[1], @srcs),") \$(LDFLAGS) -o \$@ \$(OBJS) \$(LIBS)\n";
print MAKEFILE "\t\@AO1=`mktemp _ao_XXXXXX`;\\\n";
print MAKEFILE "\tAO2=`mktemp _ao_XXXXXX`;\\\n";
print MAKEFILE "\techo \$(OBJS) | tr ' ' '\\n' | sort > \$AO1;\\\n";
print MAKEFILE "\techo *.o     | tr ' ' '\\n' | sort > \$AO2;\\\n";
print MAKEFILE "\tif [ `join -v 2 \$AO1 \$AO2 | wc -l` -gt 0 ] ; then\\\n";
print MAKEFILE "\t\techo -n \"WARNING: unused object files: \";\\\n";
print MAKEFILE "\t\tjoin -v 2 \$AO1 \$AO2 | tr '\\n' ' ';\\\n";
print MAKEFILE "\t\techo;\\\n";
print MAKEFILE "\tfi;\\\n";
print MAKEFILE "\t\$(RM) \$AO1 \$AO2\n\n";
#
# make date
#
print MAKEFILE "date: \n";
print MAKEFILE "\t\@if [ ! -f env.dat ]; then\\\n";
print MAKEFILE "\t\techo \"". $setup_line . "\" > env.dat; \\\n";
print MAKEFILE "\t\tsed -n '/Id:/p' *.h *.c *.F90 | column -t >> env.dat; \\\n";
print MAKEFILE "\t\tsed -e '/^\$\$/ d' -e \"/^\\// d\" piernik.def >> env.dat; \\\n";
print MAKEFILE "\tfi;\n\n";
#
# make version.F90
#
print MAKEFILE "version.F90: date\n";
print MAKEFILE "\t\@if [ -e version.F90 ]; then unlink version.F90; fi; \n";
print MAKEFILE "\t@( echo -e \"module version\\n   implicit none\\n   public\\n\"; \\\n";
print MAKEFILE "\twc -l env.dat | awk '{print \"   integer, parameter :: nenv = \"\$\$1\"+0\"}'; \\\n";
print MAKEFILE "\techo -e \"   character(len=128), dimension(nenv) :: env\\ncontains\\n   subroutine init_version\\n      implicit none\"; \\\n";
print MAKEFILE "\tawk '{printf(\"       env(%i) = \\\"%s\\\"\\n\",NR,\$\$0)}' env.dat; \\\n";
print MAKEFILE "\techo -e \"    end subroutine init_version\\nend module version\" ) > version.F90; \\\n";
print MAKEFILE "\techo 'generated version.F90'; \\\n\n";
#
# make clean
#
print MAKEFILE "clean:\n";
print MAKEFILE "\t\$(RM) \$(PROG) \$(OBJS) *.mod\n\n";
#
# make clean-run
#
print MAKEFILE "clean-run:\n";
print MAKEFILE "\t\$(RM) *.bck *~ *.h5 *.res *.log *.tsl *.out *.tmp core*\n\n";
#
# make clean-all
#
print MAKEFILE "clean-all:\n";
print MAKEFILE "\t\$(RM) \$(PROG) \$(OBJS) *.mod *.bck *~ *.h5 *.res *.log *.tsl *.out *.tmp core* *.f *.dbg \n\n";
#
# Make .F90 a valid suffix
#
print MAKEFILE ".SUFFIXES: \$(SUFFIXES) .F90\n\n";
#
# .F90 -> .o
#
print MAKEFILE ".F90.o:\n";
print MAKEFILE "\t\$(ECHO_FC)\n";
print MAKEFILE "ifdef USE_GNUCPP\n";
print MAKEFILE "\t\@\$(PRECOMP) \$(CPPFLAGS) \$< \$(patsubst %.F90,%\_.f90,\$<)\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c \$(patsubst %.F90,%\_.f90,\$<) && ( \$(RM) \$(patsubst %.F90,%\_.f90,\$<); \$(MV) \$(patsubst %.F90,%\_.o,\$<) \$(patsubst %.F90,%.o,\$<) )\n";
print MAKEFILE "else\n";
print MAKEFILE "\t\$(F90) \$(CPPFLAGS) \$(F90FLAGS) -c \$<\n";
print MAKEFILE "endif\n\n";
#
# .c -> .o
#
print MAKEFILE ".c.o:\n";
print MAKEFILE "\t\$(ECHO_CC)\n";
print MAKEFILE "\t\$(CC) \$(CPPFLAGS) \$(CFLAGS) -c \$<\n\n";
#
# override the built-in rule for .mod (Modula-2 source code files)
#
print MAKEFILE "%.o : %.mod\n\n";
#
# Dependency listings
#
&MakeDependsf90($ARGV[1]);
&MakeDepends("*.f *.F *.F90", '^\s*use\s+["\']([^"\']+)["\']');
&MakeDepends("*.f *.F *.F90", '^\s*#\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.c",     '^\s*#\s*include\s+["\']([^"\']+)["\']');

# system("./post","mhd",$ARGV[0]); # !found
#system("./newcompiler Makefile");
print "Setting the '" . $compiler . "' configuration.\n\n";
system("./newcompiler " . substr($compiler, 0, -3) );

# COMPILE NOW
chdir('..');
# do checks
if ( (! -r $objdir) || (! -r "$objdir/".$param_file) || (! -r "$objdir/Makefile") ) {
   print "problem with directory \'$objdir\'\n";
   exit $FAILED;
}

if ($moo{nocompile}) {
  print "Compilation skipped on request. You may want to run 'make -C $objdir' to compile the Piernik code.\n\n";
  @cir = grep(/Circular/ig, `LC_ALL=C make -C $objdir CHECK_MAGIC=yes piernik 2>&1`);
  if( @cir ) {
     print @cir;
     exit $FAILED
  }
  exit $OK;
}

# compile
system("make -C $objdir");

if ( ! -x "$objdir/piernik" ) {
   print "\nit appears that \'make\' crashed\n";
   print "cannot continue\n";
   exit $FAILED;
}

my $problemname = $ARGV[0];
if ( -d 'run' ) {
    print "\ndirectory \'run\' exists but not used anymore\n";
    print "point your scripts to \'runs/" . $problemname . "\' instead\n";
}
if ( ! -d 'runs' ) {
    mkdir('runs');
}
my $problempath = "runs/" . $problemname;
if ($moo{obj}) { $problempath .= "_$moo{obj}"; }
if ( ! -d $problempath ) {
    mkdir($problempath);
} else {
  if ($moo{obj}) {
    unlink "$problempath/piernik", "$problempath/piernik.def";
    rename "$problempath/problem.par", "$problempath/problem.par.old";
  } else {
    rmtree([$problempath]);
    mkpath([$problempath]);
  }
}

my $runOK;
if (copy("$objdir/$param_file", "$problempath/problem.par") && \
    copy("$objdir/piernik.def", $problempath) ) {
  if ($moo{linkexe}) { $runOK = symlink("../../$objdir/piernik", "$problempath/piernik"); }
  else               { $runOK = copy("$objdir/piernik", $problempath) && chmod(0755, "$problempath/piernik"); }
}
if ($runOK) { print "\n $problemname ready in $problempath\n"; }
else {
   print "$!\n";
   exit $FAILED;
}

exit $OK;

#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "FC"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CC"; last CASE; };
         $compiler = "F90";
         }
      }
   else {
      CASE: {
         grep(/\.F90$/, @srcs)   && do { $compiler = "F90"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "FC";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
         $compiler = "???";
         }
      }
   $compiler;
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, $1);
         }
      if (defined @incs) {
         $file =~ s/\.[^.]+$/.o/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @incs);
         print MAKEFILE "\n";
         undef @incs;
         }
      }
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.F90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.F90$/.o/;
         }
      }
   $filename{"version"} = "version.o";
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.F90>) {
      `cpp -w $cppflags $file > $file"_"`;
      open(FILE, $file."_");
      while (<FILE>) {
#         /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      if (defined @incs || defined @modules) {
         ($objfile = $file) =~ s/\.F90$/.o/;
         print MAKEFILE "$objfile: $file ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies, $filename{$module});
            }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         #
         # Cray F90 compiler
         #
         if ($compiler eq "cray") {
            print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c ";
            foreach $depend (@dependencies) {
               push(@modules, "-p", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         #
         # ParaSoft F90 compiler
         #
         if ($compiler eq "parasoft") {
            print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c ";
            foreach $depend (@dependencies) {
               $depend =~ s/\.o$/.F90/;
               push(@modules, "-module", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         }
      unlink $file."_";
      }
   }
