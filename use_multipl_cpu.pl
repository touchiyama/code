#!/home/touchiyama/.linuxbrew/bin/perl5.30.0

use lib ('/home/touchiyama/perl5/lib/perl5/Parallel/');
use Parallel::ForkManager;

$file=$ARGV[0];
$N_cpu=$ARGV[1];

open(IN,$file);
while(my $line=<IN>){
    chomp($line);
    push(@run,$line);
}
close(IN);

$pm=new Parallel::ForkManager($N_cpu);
foreach my $run (@run){
    my $pid=$pm -> start and next;
    system "$run";
    $pm->finish;
}
