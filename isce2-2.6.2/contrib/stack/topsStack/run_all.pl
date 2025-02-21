#!/usr/bin/perl
use FileHandle;
$log = "";

if (($#ARGV + 1) < 2){die <<EOS ;}
*** run_all.pl: Run a single command interating over arguments stored in a list, each line of the list provides a set of arguments for the command ***
*** Copyright 2019, Gamma Remote Sensing, v1.3 22-Mar-2019 clw/cm ***

usage: $0 <list> <command> [log] 
    list      (input) Multi-column list of arguments for the command
              Note: Each line of the list contains a set of arguments \$1 \$2... that are inserted into the template and then executed
    command   Command template, must be entered in single quotes, placeholders for command arguments supplied from the list are specified as \$1, \$2 ...
              Example:  'cp -r \$1 \$2'  where \$1 is the source file and \$2 is the destination file or directory.
              In this example each row of the list supplies 2 new values for the source and destination files or directories
    log       (output) optional log file that captures all screen output, both stdout and stderr are directed to the log file

    Example: run_all.pl mylist 'cp -r \$1 \$2' cp.log
    
EOS

$ltab  = $ARGV[0];
$tplt  = $ARGV[1];
open(LTAB, "<$ltab") or die "\nERROR $0: list of arguments does not exist: $ltab\n\n";

$time = localtime;
print "start: $time\n";
print "command: \'$tplt\'\n";
print "command list: $ltab\n";

if($#ARGV >= 2){
  $log = $ARGV[2];
  print "log file: $log\n\n";
  open(LOG,">$log") or die "ERROR $0: cannot open log file: $log\n";
  print LOG "start: $time\n";
  print LOG "log file: $log\n";
  print LOG "command: \'$tplt\'\n";
  print LOG "command list: $ltab\n\n";  
}
else {print "\n";}

$i = 1;
LINE: while (<LTAB>) {#read lines of processing list file
  chomp $_;		#remove new line from record
  next LINE if /^$/; 	#skip blank lines in command line argument list
  next LINE if /^#/; 	#skip comments in command line argument list
  @fields = split;
  $nv = scalar @fields;
  $cmd = $tplt;

  for ($j=1; $j <= $nv; $j++) {
    $srch ="\$$j";
    $cmd =~ s/\Q$srch\E/$fields[$j-1]/g;  #regex escape $ in the search string with \Q \E, replace all instances of $1, $2 ... with /g option
#    print "search: $srch  replace: $fields[$j-1]  cmd: $cmd\n";
  }

  print "iteration: $i\n$cmd\n";
  $result = `$cmd 2>&1`;  	#capture both STDOUT and STDERR
  print "$result";
  if ($log ne ""){print LOG "$result\n";}
  $i++;
}

$time = localtime;
print "$0 completed: $time\n";
if ($log ne ""){
  print LOG "$0 completed: $time\n";
  close LOG;
}
exit 0;
