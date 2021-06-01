#!/usr/bin/perl

use lib ('/home00/kohanada/bin');
use ForkManager;

$file=$ARGV[0];
$N_cpu=$ARGV[1];
open(in, $file);
while(<in>){
        chomp;
        push(@run, $_);
}
close in;

$pm=new Parallel::ForkManager($N_cpu);
for(my $i=0; $i<@run; $i++){
        my $pid=$pm -> start and next;
        `$run[$i]`;
        $pm->finish;
}

