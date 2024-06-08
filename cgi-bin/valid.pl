#!/usr/bin/perl

############################### Header Information ##############################
require 'cgi.perl';
use CGI;;
$query = new CGI;
&ReadParse;
print &PrintHeader;

################################ Reads Inoput Data ##############################
$atom = $query->param('atom');
$file = $query->param('file');
$svm_th = $query->param('svm_th');

#################Validation Of Input Sequence Data (file upload) ###################################
if($file ne '' && $atom eq '')
{
    $file=~m/^.*(\\|\/)(.*)/; 
    while(<$file>) 
    {
	$seqfi .= $_;
    }
}
elsif($atom ne '' && $file eq ''){

    $seqfi="$atom";
}

##############ACTUAL PROCESS BEGINS FROM HERE#######################
$infut_file = "/webservers/cgi-bin/subnucpred";
$ran= int(rand 10000);
$dir = "/webservers/cgidocs/mkumar/temp/Ravindra/SubNucPred/subnuc$ran";
#$dir = "/webservers/cgidocs/mkumar/temp/Ravindra/SubNucPred/subnuc$ran";
system "mkdir $dir";
system "chmod 777 $dir";
$nam = 'input.'.'fasta';
open(FP1,">$dir/input_meta.fasta");
print FP1 "$seqfi\n";
#print "$seqfi\n";
close FP1;

system "/usr/bin/tr -d '\r' <$dir/input_meta.fasta >$dir/input_fi.fasta";
system "/usr/bin/perl $infut_file/fasta.pl $dir/input_fi.fasta |/usr/bin/head -50 >$dir/input.fasta";
system "/bin/grep -c '>' $dir/input.fasta >$dir/total_seq";
system "/bin/grep '>' $dir/input.fasta |/usr/bin/cut -d '|' -f3 |/usr/bin/cut -d ' ' -f1 >$dir/protein_id";
system "/usr/local/bin/hmmscan --domtblout $dir/hmm_out -E 1e-5 $infut_file/Pfam/Pfam-A.hmm $dir/input.fasta >/dev/null";
system "/usr/bin/perl $infut_file/hmm.pl $dir/hmm_out |/bin/grep -v '#' >$dir/hmm_out_domain";
system "/usr/bin/perl $infut_file/evalue.pl $dir/hmm_out_domain |/usr/bin/cut -d ' ' -f1,2,3 |/usr/bin/tr '|' '#' >$dir/hmm_out_domain_evalue";
system "/usr/bin/perl $infut_file/id_compare.pl $dir/hmm_out_domain_evalue $dir/protein_id |/usr/bin/cut -d ' ' -f1,2 |/usr/bin/sort -u >$dir/domain_id";
system "/usr/bin/perl $infut_file/aaseqformat.pl $dir/input.fasta |/bin/sed -e 's/^/+1/' >$dir/aacomp";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_cent $dir/svm_score_cent >/dev/null";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_chromo $dir/svm_score_chromo >/dev/null";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_nspeckle $dir/svm_score_nspeckle >/dev/null";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_nucleolus $dir/svm_score_nucleolus >/dev/null";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_other $dir/svm_score_other >/dev/null";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_nenvelop $dir/svm_score_nenvelop >/dev/null";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_nmatrix $dir/svm_score_nmatrix >/dev/null";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_npc $dir/svm_score_npc >/dev/null";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_nplasm $dir/svm_score_nplasm >/dev/null";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_pml $dir/svm_score_pml >/dev/null";
system "/usr/local/bin/svm_classify $dir/aacomp $infut_file/Models/model_tel $dir/svm_score_tel >/dev/null";
system "/usr/bin/paste $dir/protein_id $dir/aacomp $dir/svm_score_cent $dir/svm_score_chromo $dir/svm_score_nspeckle $dir/svm_score_nucleolus $dir/svm_score_other $dir/svm_score_nenvelop $dir/svm_score_nmatrix $dir/svm_score_npc $dir/svm_score_nplasm $dir/svm_score_pml $dir/svm_score_tel |/usr/bin/tr '\t' '#' >$dir/final";
system "/usr/bin/perl $infut_file/domain.pl $dir/final $dir/domain_id >$dir/final_pred";
$total_seq=`head -1 $dir/total_seq`;chomp($total_seq);

open(FINAL_PRED,"$dir/final_pred") or die "$!";
@array=<FINAL_PRED>;
close FINAL_PRED;

for($q=0;$q<=$#array;$q++)
{
    @dom=split(/\#/,$array[$q]);
    $a=0;$b=0;$c=0;$d=0;$e=0;$f=0;$g=0;$h=0;$i=0;$j=0;
    open(UNIQUE_DOMAIN,"$infut_file/nuclear_unique_domain") or die "$!";
    while($domainfile=<UNIQUE_DOMAIN>)
    {
	chomp($domainfile);
	@domain=split(/\#/,$domainfile);
	for($x=13;$x<=$#dom;$x++)
	{
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "CENTROMERE"))
	    {
		$a=100;
	    }
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "CHROMOSOME"))
	    {
		$b=100;
	    }
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "NSPECKLE"))
	    {
		$c=100;
	    }
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "NUCLEOLUS"))
	    {
		$d=100;
	    }
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "NENVELOP"))
	    {
		$e=100;
	    }
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "NMATRIX"))
	    {
		$f=100;
	    }
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "NPC"))
	    {
		$g=100;
	    }
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "NPLASM"))
	    {
		$h=100;
	    }
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "PML"))
	    {
		$i=100;
	    }
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "TELOMERE"))
	    {
		$j=100;
	    }
	}
    }
    close UNIQUE_DOMAIN;
    if($a==100)
    {
	open(RESULT,">>$dir/Hybrid_result") or die "$!";
	print RESULT "$dom[0]#CENTROMERE#Domain_based_prediction\n";
	close RESULT;
    }
    if($b==100)
    {
	open(RESULT,">>$dir/Hybrid_result") or die "$!";
	print RESULT "$dom[0]#CHROMOSOME#Domain_based_prediction\n";
	close RESULT;
    }
    if($c==100)
    {
	open(RESULT,">>$dir/Hybrid_result") or die "$!";
	print RESULT "$dom[0]#NUCLEAR_SPECKLE#Domain_based_prediction\n";
	close RESULT;
    }
    if($d==100)
    {
	open(RESULT,">>$dir/Hybrid_result") or die "$!";
	print RESULT "$dom[0]#NUCLEOLUS#Domain_based_prediction\n";
	close RESULT;
    }
    if($e==100)
    {
	open(RESULT,">>$dir/Hybrid_result") or die "$!";
	print RESULT "$dom[0]#NUCLEAR_ENVELOPE#Domain_based_prediction\n";
	close RESULT;
    }
    if($f==100)
    {
	open(RESULT,">>$dir/Hybrid_result") or die "$!";
	print RESULT "$dom[0]#NUCLEAR_MATRIX#Domain_based_prediction\n";
	close RESULT;
    }
    if($g==100)
    {
	open(RESULT,">>$dir/Hybrid_result") or die "$!";
	print RESULT "$dom[0]#NUCLEAR_PORE_COMPLEX#Domain_based_prediction\n";
	close RESULT;
    }
    if($h==100)
    {
	open(RESULT,">>$dir/Hybrid_result") or die "$!";
	print RESULT "$dom[0]#NUCLEOPLASM#Domain_based_prediction\n";
	close RESULT;
    }
    if($i==100)
    {
	open(RESULT,">>$dir/Hybrid_result") or die "$!";
	print RESULT "$dom[0]#PML#Domain_based_prediction\n";
	close RESULT;
    }
    if($j==100)
    {
	open(RESULT,">>$dir/Hybrid_result") or die "$!";
	print RESULT "$dom[0]#TELOMERE#Domain_based_prediction\n";
	close RESULT;
    }
    
    if(($a==0)&&($b==0)&&($c==0)&&($d==0)&&($e==0)&&($f==0)&&($g==0)&&($h==0)&&($i==0)&&($j==0))
    {
	open(ID,">>$dir/Hybrid_result") or die "$!";
	print ID "$dom[0]#CENTROMERE#$dom[2]#CHROMOSOME#$dom[3]#NUCLEAR_SPECKLE#$dom[4]#NUCLEOLUS#$dom[5]#OTHER#$dom[6]#NUCLEAR_ENVELOPE#$dom[7]#NUCLEAR_MATRIX#$dom[8]#NUCLEAR_PORE_COMPLEX#$dom[9]#NUCLEOPLASM#$dom[10]#PML body#$dom[11]#TELOMERE#$dom[12]\n";
	close ID;
    }
}
open(NE,"$dir/Hybrid_result") or die "$!";
while($tim=<NE>)
{
    chomp($tim);
    if($tim =~ m/Domain_based_prediction$/)	 
    {
	open(YY,">>$dir/result_both") or die "$!";
	print YY "$tim\n";
	close YY;
    }
    if($tim!~ m/Domain_based_prediction$/)
    {
	@km=split(/\#/,$tim);
	open(YY,">>$dir/result_both") or die "$!";
	print YY "$km[0]#";
	close YY;
	if($km[2] >= $svm_th)
	{
	    open(YY,">>$dir/result_both") or die "$!";
	    print YY "$km[1]#$km[2]#";
	    close YY;
	}
	if($km[4] >= $svm_th)
	{
	    open(YY,">>$dir/result_both") or die "$!";
	    print YY "$km[3]#$km[4]#";
	    close YY;
	}
	if($km[6] >= $svm_th)
	{
	    open(YY,">>$dir/result_both") or die "$!";
	    print YY "$km[5]#$km[6]#";
	    close YY;
	}
	if($km[8] >= $svm_th)
	{
	    open(YY,">>$dir/result_both") or die "$!";
	    print YY "$km[7]#$km[8]#";
	    close YY;
	}
	if(($km[10] >= $svm_th)&&($km[12] >= $svm_th))
	{
	    open(YY,">>$dir/result_both") or die "$!";
	    print YY "$km[11]#$km[12]#";
	    close YY;
	}
	if(($km[10] >= $svm_th)&&($km[14] >= $svm_th))
	{
	    open(YY,">>$dir/result_both") or die "$!";
	    print YY "$km[13]#$km[14]#";
	    close YY;
	}
	if(($km[10] >= $svm_th)&&($km[16] >= $svm_th))
	{
	    open(YY,">>$dir/result_both") or die "$!";
	    print YY "$km[15]#$km[16]#";
	    close YY;
	}
	if(($km[10] >= $svm_th)&&($km[18] >= $svm_th))
	{
	    open(YY,">>$dir/result_both") or die "$!";
	    print YY "$km[17]#$km[18]#";
	    close YY;
	}
	if(($km[10] >= $svm_th)&&($km[20] >= $svm_th))
	{
	    open(YY,">>$dir/result_both") or die "$!";
	    print YY "$km[19]#$km[20]#";
	    close YY;
	}
	if(($km[10] >= $svm_th)&&($km[22] >= $svm_th))
	{
	    open(YY,">>$dir/result_both") or die "$!";
	    print YY "$km[21]#$km[22]#";
	    close YY;
	}
	open(YY,">>$dir/result_both") or die "$!";
	print YY "\n";
	close YY;
    }
}

print  "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n";
print  "<html><HEAD>\n";
print  "<TITLE>SubNucPred::Prediction Result</TITLE>\n";
print  "<META NAME=\"description\" CONTENT=\"SubNucPred, University of Delhi South Campus, INDIA\">\n";
print  "</HEAD><body bgcolor=\"\#FFFFE0\">\n";
print  "<h2 ALIGN = \"CENTER\"> SubNucPred Prediction Result</h2>\n";
print  "<HR ALIGN =\"CENTER\"> </HR>\n";
if($total_seq == 0)
{
    print  "<p align=\"center\"><font size=4 color=black><b>Please submit your sequence in fasta format</b></p>";
}

else
{
    print  "<p align=\"center\"><font size=4 color=black><b>The submitted protein/proteins belongs to <font color='red'></p>";
    print "<table border='1' width='544' align='center'>";
    print "<tr><th>Protein ID</th><th>Location</th><th>Score</th></tr>";

    open(PREDICTION,"$dir/result_both") or die "$!";
    while($pre=<PREDICTION>)
    {
	chomp($pre);
	@pred=split(/\#/,$pre);
	$no=0;
	for($kk=0;$kk<$#pred;$kk++)
	{
	    $no++;
	}
	$mo=$no/2;
	if($pre =~ m/Domain_based_prediction$/)
	{
	    print "<tr align='center' ><td>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td></tr>";
	}
	if($pre!~ m/Domain_based_prediction$/)
	{
	    if($no==0)
	    {
		print "<tr align='center' ><td rowspan='$mo'>$pred[0]</td><td>No_Prediction</td><td></td></tr>";
	    }
	    if($no==2)
	    {
		print "<tr align='center' ><td rowspan='$mo'>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td></tr>";
	    }
	    if($no==4)
	    {
		print "<tr align='center' ><tr align='center'><td rowspan='$mo'>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td><tr align='center'><td>$pred[3]</td><td>$pred[4]</td></tr></tr></tr>";
	    }
	    if($no==6)
	    {
		print "<tr align='center' ><tr align='center'><td rowspan='$mo'>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td><tr align='center'><td>$pred[3]</td><td>$pred[4]</td><tr align='center'><td>$pred[5]</td><td>$pred[6]</td></tr></tr></tr></tr>";
	    }
	    if($no==8)
	    {
		print "<tr align='center' ><tr align='center'><td rowspan='$mo'>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td><tr align='center'><td>$pred[3]</td><td>$pred[4]</td><tr align='center'><td>$pred[5]</td><td>$pred[6]</td><tr align='center'><td>$pred[7]</td><td>$pred[8]</td></tr></tr></tr></tr></tr>";
	    }
	    if($no==10)
	    {
		print "<tr align='center' ><tr align='center'><td rowspan='$mo'>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td><tr align='center'><td>$pred[3]</td><td>$pred[4]</td><tr align='center'><td>$pred[5]</td><td>$pred[6]</td><tr align='center'><td>$pred[7]</td><td>$pred[8]</td><tr align='center'><td>$pred[9]</td><td>$pred[10]</td></tr></tr></tr></tr></tr></tr>";
	    }
	    if($no==12)
	    {
		print "<tr align='center' ><tr align='center'><td rowspan='$mo'>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td><tr align='center'><td>$pred[3]</td><td>$pred[4]</td><tr align='center'><td>$pred[5]</td><td>$pred[6]</td><tr align='center'><td>$pred[7]</td><td>$pred[8]</td><tr align='center'><td>$pred[9]</td><td>$pred[10]</td><tr align='center'><td>$pred[11]</td><td>$pred[12]</td></tr></tr></tr></tr></tr></tr></tr>";
	    }
	    if($no==14)
	    {
		print "<tr align='center' ><tr align='center'><td rowspan='$mo'>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td><tr align='center'><td>$pred[3]</td><td>$pred[4]</td><tr align='center'><td>$pred[5]</td><td>$pred[6]</td><tr align='center'><td>$pred[7]</td><td>$pred[8]</td><tr align='center'><td>$pred[9]</td><td>$pred[10]</td><tr align='center'><td>$pred[11]</td><td>$pred[12]</td><tr align='center'><td>$pred[13]</td><td>$pred[14]</td></tr></tr></tr></tr></tr></tr></tr></tr>";
	    }
	    if($no==16)
	    {
		print "<tr align='center' ><tr align='center'><td rowspan='$mo'>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td><tr align='center'><td>$pred[3]</td><td>$pred[4]</td><tr align='center'><td>$pred[5]</td><td>$pred[6]</td><tr align='center'><td>$pred[7]</td><td>$pred[8]</td><tr align='center'><td>$pred[9]</td><td>$pred[10]</td><tr align='center'><td>$pred[11]</td><td>$pred[12]</td><tr align='center'><td>$pred[13]</td><td>$pred[14]</td><tr align='center'><td>$pred[15]</td><td>$pred[16]</td></tr></tr></tr></tr></tr></tr></tr></tr></tr>";
	    }
	    if($no==18)
	    {
		print "<tr align='center' ><tr align='center'><td rowspan='$mo'>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td><tr align='center'><td>$pred[3]</td><td>$pred[4]</td><tr align='center'><td>$pred[5]</td><td>$pred[6]</td><tr align='center'><td>$pred[7]</td><td>$pred[8]</td><tr align='center'><td>$pred[9]</td><td>$pred[10]</td><tr align='center'><td>$pred[11]</td><td>$pred[12]</td><tr align='center'><td>$pred[13]</td><td>$pred[14]</td><tr align='center'><td>$pred[15]</td><td>$pred[16]</td><tr align='center'><td>$pred[17]</td><td>$pred[18]</td></tr></tr></tr></tr></tr></tr></tr></tr></tr></tr>";
	    }
	    if($no==20)
	    {
		print "<tr align='center' ><tr align='center'><td rowspan='$mo'>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td><tr align='center'><td>$pred[3]</td><td>$pred[4]</td><tr align='center'><td>$pred[5]</td><td>$pred[6]</td><tr align='center'><td>$pred[7]</td><td>$pred[8]</td><tr align='center'><td>$pred[9]</td><td>$pred[10]</td><tr align='center'><td>$pred[11]</td><td>$pred[12]</td><tr align='center'><td>$pred[13]</td><td>$pred[14]</td><tr align='center'><td>$pred[15]</td><td>$pred[16]</td><tr align='center'><td>$pred[17]</td><td>$pred[18]</td><tr align='center'><td>$pred[19]</td><td>$pred[20]</td></tr></tr></tr></tr></tr></tr></tr></tr></tr></tr></tr>";
	    }
	}
    }
}
print "</table>";
print "</font></b></font></p>\n";
print  "<p align=\"center\"><font size=3 color=black><b>Thanks for using SubNucPred Prediction Server</b></font></p>\n";
print  "<p align=\"center\"><font size=3 color=black><b>If you have any problem or suggestions please contact <a href='mailto:manish@south.du.ac.in'>Dr. Manish Kumar</a></b></font>. Please mention your job number in any communication.</p></br>\n";

if($total_seq >= 1)
{
    print  "<p ALIGN=\"CENTER\"><b>Your job number is <font color=\"red\">$ran</b></font></p>\n";
}
print  "</body>\n";
print  "</html>\n";
