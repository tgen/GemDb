#!/usr/bin/perl 
package updateGeneLocat;
use initMongo;

sub loadTo {
    my $fileInsert = initMongo->addToCollections(@_);
    my $function   = shift;
    my $filePath   = shift;
    my $collection = shift;
    my $hostName   = shift;
    my $key        = shift;
    print
"updateGeneLocat: key: $key filepath:$filePath collection:$collection hostName: $hostName\n";
    my %done = ();
    my ( $conn, $db, $GENE ) =initMongo->new( 'AnnotateBy', $collection, $hostName );
    my ( $connMarkers, $dbMarkers, $MARKERS ) = initMongo->new( 'markers', 'tumor', $hostName );
    my ( $connGermMarkers, $GermDbMarkers, $GERMMARKERS ) = initMongo->new( 'markers', 'germline', $hostName );    
    open( DATA, "$filePath" ) or die "Can't open $filePath\n";
    LOOP: while (<DATA>) {
        if (/^#/) { next LOOP }
        chomp;
        $line = $_;
        @temp = split( /\t/, $line );
        if ( $temp[2] eq "gene" ) {
            if ( $temp[8] =~ /gene_biotype \"protein_coding"/ ) {
                my $geneName = "none";
                if ( $temp[8] =~ /gene_name \"(.*?)\"/ ) {
                    $geneName = $1;
                    unless ($temp[0]=~/^G/) {next LOOP;}
                    $chr      = 'chr' . $temp[0];
                    $startPos = int( $temp[3] );
                    $endPos   = int( $temp[4] );
                    $dir      = $temp[6];
                    print "geneName:$geneName - $chr $startPos $endPos $dir\n";
                    $GENE->update({"gene" => $geneName},{'$set'=>{'startPos'=>$startPos,'chr'=>$chr,'endPos'=>$endPos,'dir'=>$dir}},{'multiple'=>1,'safe'=>1});
                    $MARKERS->update({"gene" => $geneName},{'$set'=>{'geneInfo.startPos'=>$startPos,'geneInfo.chr'=>$chr,'geneInfo.endPos'=>$endPos,'geneInfo.dir'=>$dir}},{'multiple'=>1,'safe'=>1});
                    $GERMMARKERS->update({"gene" => $geneName},{'$set'=>{'geneInfo.startPos'=>$startPos,'geneInfo.chr'=>$chr,'geneInfo.endPos'=>$endPos,'geneInfo.dir'=>$dir}},{'multiple'=>1,'safe'=>1});
                }
            }
        }
    }
}
1;
