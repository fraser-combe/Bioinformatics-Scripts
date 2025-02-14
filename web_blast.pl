#!/usr/bin/perl

use strict;
use warnings;
use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);

$| = 1;

my $ua = LWP::UserAgent->new;
$ua->agent("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/89.0.4389.114 Safari/537.36");

my $fasta_file = shift @ARGV or die "Usage: $0 <fasta_file>\n";

# Read and encode the sequences from the FASTA file
open(my $fh, '<', $fasta_file) or die "Could not open file '$fasta_file' $!";
my $encoded_query = "";
while (my $line = <$fh>) {
    $encoded_query .= uri_escape($line);
}
close($fh);

# Build the request
my $args = "CMD=Put&PROGRAM=blastn&DATABASE=nt&QUERY=" . $encoded_query;

my $req = HTTP::Request->new(POST => 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi');
$req->content_type('application/x-www-form-urlencoded');
$req->content($args);

# Get the response
my $response = $ua->request($req);

# Log the response content to check for issues
open(my $log_fh, '>', 'blast_log.txt') or die "Could not open file 'blast_log.txt' $!";
print $log_fh "Response content:\n" . $response->content . "\n";
close($log_fh);

# Parse out the request id
$response->content =~ /^    RID = (.*$)/m;
my $rid = $1;

unless ($rid) {
    die "Failed to retrieve RID from BLAST response. Check blast_log.txt for details.";
}

# Parse out the estimated time to completion
$response->content =~ /^    RTOE = (.*$)/m;
my $rtoe = $1;

# Wait for search to complete
sleep $rtoe;

# Poll for results
my $timeout = 300; # 5 minutes timeout
my $elapsed = 0;
while (true) {
    sleep 5;
    $elapsed += 5;

    last if $elapsed > $timeout;

    $req = HTTP::Request->new(GET => "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid");
    $response = $ua->request($req);

    if ($response->content =~ /\s+Status=WAITING/m) {
        next;
    }

    if ($response->content =~ /\s+Status=FAILED/m) {
        die "Search $rid failed; please report to blast-help\@ncbi.nlm.nih.gov.\n";
    }

    if ($response->content =~ /\s+Status=UNKNOWN/m) {
        die "Search $rid expired.\n";
    }

    if ($response->content =~ /\s+Status=READY/m) {
        if ($response->content =~ /\s+ThereAreHits=yes/m) {
            last;
        } else {
            die "No hits found.\n";
        }
    }
}

# Retrieve and display results
$req = HTTP::Request->new(GET => "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=$rid");
$response = $ua->request($req);

print $response->content;
exit 0;
