#!/usr/bin/env perl

# Code to pull out the dependency tree for .tex files, knows about pdf figures

$docpath=$ENV{"PWD"};
@Current = @ARGV;

until ($end) {
    push(@Depending,&dependants(@Current));
    push(@graphics,&graphics(@Current));
    
    local(%MARKER);
    grep($MARKER{$_}++,@Current);
    @clash=grep($MARKER{$_},@Depending);
    if (@clash) {die "dependence heading immediately up the document tree!\n@clash\n"}
    undef @clash;

    local(%MARKER);
    grep($MARKER{$_}++,@Depending);
    @clash=grep($MARKER{$_},@Previous);
    if (@clash) {die "dependence heading up the document tree!\n@clash\n"}
    
    push(@Previous,@Current);
    undef @Current;

    if (@Depending) {
	push(@Current,@Depending);
	undef @Depending;
    } else {
	$end=1;
    }
}

print "MANUALCOMPONENTS = ";
foreach (@Previous) {print "$_ "}
foreach (@graphics) {print "$_ "}
print "\n";

sub dependants {
    my (@files) = @_;
    my (@texfiles);
    foreach $file (@files) {
	open(FILE,"$docpath/$file") || die "Can't open $docpath/$file\n";
	my @contents = <FILE>;
	close(FILE);
	foreach (@contents) {
	    if (/^\s*\\input\{(.+)\.?t?e?x?\}/ && !/^\s*\%/) {
		my $texfile = $1;
		$texfile =~ s/\.tex$//;
		push(@texfiles,"$texfile\.tex");
	    }
	    if (/^\s*\\include\{(.+)\.?t?e?x?\}/ && !/^\s*\%/) {
		my $texfile = $1;
		$texfile =~ s/\.tex$//;
		push(@texfiles,"$texfile\.tex");
	    }
	}
    }
    @texfiles;
}

sub graphics {
    my (@files) = @_;
    my (@pdffiles);
    foreach $file (@files) {
        open(FILE,"$docpath/$file") || die "Can't open $docpath/$file\n";
        my @contents = <FILE>;
        close(FILE);
        foreach (@contents) {
            if (/^\s*\\includegraphics.*\{(.+)\}/ && !/^\s*\%/) {
                push(@pdffiles,$1);
            }
        }
    }
    @pdffiles;
}
