my @fiels;
my $line;
my $score=100;
my $id="";

$line=<>;
@fields=split(/\t/,$line);
$id=$fields[0];
$score=$fields[10];

while(<>) {
        @fields=split(/\t/,$_);
        if ($fields[0] ne $id) {
                $id=$fields[0];
                print $line;
                $score=100;
        }
        if ($score>$fields[10]) {
                $score=$fields[10];
                $line=$_;
        }
}



