##!/usr/bin/perl

#Author: huangjingfeng@genomics.org.cn
#Date: Thu Oct 28 15:59:46 CST 2016
#use strict;
#use warnings;
#use Data::Dumper;
use FindBin '$Bin';
use Getopt::Long;
use threads;

my ($diffdir,$explist,$colgeneid,$udcol,$colexp,$outdir,$config,$symbol,$help);
GetOptions(
        "diffdir=s"   => \$diffdir,
	"explist=s"   => \$explist,
	"udcol:s"  => \$udcol,
	"config=s"    => \$config,
	"symbol:s"    => \$symbol,
	"colgeneid:i" => \$colgeneid,
	"colexp:i"    => \$colexp,
	"outdir:s"    => \$outdir,
	"Rpath:s" => \$Rscript,
	"convert:s" => \$Convert,
	"help|?" => \$help
);

$outdir  ||= "./";
#$symbol  ||= "/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2016a/Database/hg19/human.gene2symbol.txt";
$colgeneid ||= 1;
$colexp ||= 5;
$udcol ||= 7;
#$Rscript    ||= "/ifswh1/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/software/R-3.1.1/bin/Rscript";
$Rscript    ||= "$Bin/../software/Rscript";
$Convert    ||= "$Bin/../software/convert";

if (!defined $diffdir || !defined $explist || !defined $config ||defined $help) {

        die << "USAGE";
description: cluster
usage: perl $0 [options]
options:
        -diffdir  <path>*           input result dir,containing file of diff result, e.g.,*.GeneDiffExpFilter.xls
	-explist  <file>*	    input GeneExp list,containing path of gene expression, e.g.,***.fpkm.xls
	-config	  <file>*	    config for DiffExp, e.g.,GeneDiffExp_Allin.conf
	-udcol   		    column of up/down
	-symbol			    (if you select T) input symbol2gene. default:"/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2016a/Database/hg19/human.gene2symbol.txt"
	-colgeneid                  column of geneid in GeneExp file. default: 1
	-colexp 		    column of exp in GeneExp file. default: 5
	-Rpath			    path of Rscript, default: /opt/blc/genome/bin/Rscript
	-convert		    path of convert, default: /usr/bin/convert
	-help|?                     information for help

e.g.:
        perl $0 -diffdir /ifs4/BC_RD/USER/huangjingfeng/RNApipeline_tset/DEGs/example/result/DEseq2 -explist /ifs4/BC_RD/USER/huangjingfeng/RNApipeline_tset/DEGs/example/list/GeneExp_RSEM.downstream.list -config /ifs4/BC_RD/USER/huangjingfeng/RNApipeline_tset/DEGs/example/config/GeneDiffExp_Allin.conf  -symbol /ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2016a/Database/hg19/human.gene2symbol.txt -udcol 7 -colgeneid 1 -colexp 5 -outdir .
USAGE
}

my @diffs = glob("$diffdir/*.GeneDiffExpFilter.xls");
#die "Error: can't find file of DEGs, exit...\n" unless(@diffs);
exit unless(@diffs);

my (%exppath);
open LIST,$explist or die $!;
while (<LIST>) {
        next if (/^\s*$/ || /^\s*#/);
        chomp; my @t = split /\s+/;
        $exppath{$t[0]} = $t[1]; # exp_file   depends
}
close LIST;

my (@expline,%sample_id_exp,$expfile);
foreach my $sample (keys %exppath)
{	
	open EXP,$exppath{$sample} or die $!;
	<EXP>;
	while(<EXP>) {
	next if (/^\s*$/ || /^\s*#/);
	chomp;
	@expline = split /\t/,$_;
	my $gene_id = $expline[0];
	my $gene_exp = $expline[4];
	$sample_id_exp{$sample}{$gene_id}=$gene_exp;
	}
	close EXP;
}

my (@groupline,$method1,$dist_vs1,$method2,$dist_vs2,$dist_group,%sample_of_menthod1,%vs_of_menthod1,%vs_of_menthod2);
open CONF,$config or die $!;
$/="\n\n";
while(<CONF>)
{
	chomp;
	@groupline = split /\n/,$_;
	if (/_Group/)
	{
		if ($groupline[2]=~/(.+)_VS = (.+)/){$method1=$1;$dist_vs1=$2;}
                my @vs_element1= split /,/,$dist_vs1;
                foreach my $vs1 (@vs_element1)
                {
                        $vs_of_menthod1{$method1}{$vs1}=1;
                }
		if ($groupline[1]=~/_Group = (.+)/){$dist_group=$1;}
		my @group_element1= split /\s+/,$dist_group;
		foreach my $group_element2 (@group_element1)
		{	
			my @group2sample =split /:/,$group_element2;
			my @sample_of_group = split /,/,$group2sample[1];
			foreach my $sample1 (@sample_of_group)
			{
				$sample_of_menthod1{$method1}{$group2sample[0]}{$sample1}=1;
			}
	
		}
	}
	elsif(/_VS/) 
	{
		if ($groupline[1]=~/(.+)_VS = (.+)/){$method2=$1;$dist_vs2=$2;}
		my @vs_element2= split /,/,$dist_vs2;
		foreach my $vs2 (@vs_element2)
		{
			$vs_of_menthod2{$method2}{$vs2}=1;
		}
	}
}


#add symbol
$/="\n";
if (defined $symbol){
open SYMBOL,$symbol or die $!;
while(<SYMBOL>){
chomp;
@symbolline = split /\t/,$_;
$gene2symbol{$symbolline[0]}=$symbolline[1];
}
}

####pheatmap for vs:
print "############pheatmap for vs:\n";
my($method,$suffix,$vs,$prefix,$difffilter,@sample,@fileline,$diffgene_id,$exp_out);
foreach $method (keys %vs_of_menthod2)
{	
	if ($diffdir =~ /$method/)
	{
		$suffix = $method."_Method";
		print "$suffix:\n";
		foreach $vs (keys %{$vs_of_menthod2{$method}})
		{		
			print "$vs\n";
			@sample = split /&/,$vs;
			$prefix = "$sample[0]-VS-$sample[1]";
			####filename for script of R
			$difffilter ="$diffdir/$prefix.$suffix.GeneDiffExpFilter.xls";
			next unless(-f $difffilter);  #licong20171106
			open OUT,">$outdir/$prefix.$suffix.prepare4pheatmap-plot";
			print OUT "symbol";
			print OUT "\tUp_Down";
			close OUT;
			open OUT,">>$outdir/$prefix.$suffix.prepare4pheatmap-plot";
			foreach my $sample (@sample){print OUT "\t$sample";}
			print OUT "\n";
			#print "$difffilter\n\n";
			open IN,$difffilter or die $!;
			$/="\n";
			<IN>;
			while(<IN>)
			{
				chomp;
				@fileline = split /\t/,$_;
		 		$diffgene_id = $fileline[0];
				$up_down = $fileline[$udcol -1];
				if (exists $gene2symbol{$diffgene_id} and $gene2symbol{$diffgene_id} ne "-"){print OUT "$gene2symbol{$diffgene_id}";}
				elsif ($gene2symbol{$diffgene_id} eq "-"){print OUT "$diffgene_id";}
				else {print OUT "$diffgene_id";}
				print OUT "\t$up_down";
				foreach my $sample (@sample)
				{
					if    (!exists $sample_id_exp{$sample}{$diffgene_id}){print OUT"\t-4";}
					else  {$exp_out=log($sample_id_exp{$sample}{$diffgene_id})/log(10);print OUT"\t$exp_out";}
				}
				print OUT "\n";
			}
			close IN;
			close OUT;
			open  my $fh_rcode,">$outdir/$prefix.$suffix.pheatmap.R" or die $!;
			print $fh_rcode <<CMD;
library("pheatmap")
data <- read.table("$outdir/$prefix.$suffix.prepare4pheatmap-plot",head=T,row.names=1)
CMD
			if (!defined $symbol)
			{
print $fh_rcode <<CMD;
height <- 4
CMD
			}
			else{
print $fh_rcode <<CMD;
height <- length(row.names(data))/20
CMD
			}
print $fh_rcode <<CMD;
if (length(data) >10){
width <- 2 + length(data)*0.2
}else{
width <- 4
}
Up_Down <- data.frame(data[,1]) 
rownames(Up_Down) = paste(row.names(data))
colnames(Up_Down) = paste("Up_Down")
ann_colors = list (Up_Down=c("Up"="#E41A1C","Down"="#377EB8"))

pheatmap(
	data[,-1],
        annotation_row = Up_Down,
        annotation_names_col =F,
	annotation_names_row =F,
        annotation_colors = ann_colors, ##annotation colors
        fontsize_col=9,
	fontsize_row=2.7,
	#      border_color = "transparent",
#        color = colorRampPalette(c("#FFFFD9","yellow2","orange","orangered2","gray10"))(200), #####colors and gradient for pheatmap 
	color = colorRampPalette(c("#1086BB","#A2D4E5","#F0F4D7","#FFD694","#FD432C"))(200), #####colors and gradient for pheatmap
main="Pheatmap for $vs",
CMD
		if (!defined $symbol){
print $fh_rcode <<CMD;
	show_rownames = F,
CMD
		}
print $fh_rcode <<CMD;
	file = '$outdir/$prefix.$suffix.pheatmap-plot.pdf',
	width=width,
	height=height
)
dev.off()
CMD
		#system("export R_LIBS=$Bin/3.2:$R_LIBS && $Rscript $outdir/$prefix.$suffix.pheatmap.R 2>/dev/null");# && rm $outdir/$prefix.$suffix.pheatmap-plot.R; ###licong
		system("$Rscript $outdir/$prefix.$suffix.pheatmap.R ");
                system("$Convert -density 300 -resize 30% $outdir/$prefix.$suffix.pheatmap-plot.pdf $outdir/$prefix.$suffix.pheatmap-plot.png");
        	}
	}
}
print "done\n";


####pheatmap for group
print "##########pheatmap for group:\n";
my ($method,$suffix,$vs,@group,$prefix,$difffilter,@fileline,$diffgene_id,$exp_out);
foreach $method (keys %vs_of_menthod1)
{
	if ($diffdir =~ /$method/){
	$suffix = $method."_Method";
	print "$suffix:\n";
	foreach $vs (keys %{$vs_of_menthod1{$method}})
        {
		print "$vs\n";
		@group = split /&/,$vs;
		$prefix = "$group[0]-VS-$group[1]";
		####filename for script of R
		$difffilter ="$diffdir/$prefix.$suffix.GeneDiffExpFilter.xls";
		next unless(-f $difffilter);  #licong20171106	
		open OUT,">$outdir/$prefix.$suffix.prepare4pheatmap-plot";
                print OUT "symbol";
		print OUT "\tUp_Down";
                close OUT;
		open OUT,">>$outdir/$prefix.$suffix.prepare4pheatmap-plot";
		my @group_name;
		foreach my $group (@group)
		{
			foreach my $sample (keys %{$sample_of_menthod1{$method}{$group}})
			{
				print OUT "\t$sample";
				$group_r = join '',"\"",$group,"\"";
				push @group_name,$group_r; 
			}
		}
		my $group_array_r = join ',',@group_name;
		print OUT "\n";
		open IN,$difffilter or die $!;
		$/="\n";
		<IN>;
		while(<IN>)
		{
			chomp;
			@fileline = split /\t/,$_;
			$diffgene_id =$fileline[0];
			$up_down = $fileline[$udcol -1];
			if (exists $gene2symbol{$diffgene_id} and $gene2symbol{$diffgene_id} ne "-"){print OUT "$gene2symbol{$diffgene_id}";}
			elsif ($gene2symbol{$diffgene_id} eq "-"){print OUT "$diffgene_id";}
			else {print OUT "$diffgene_id";}
			print OUT "\t$up_down";
			foreach my $group (@group)
			{
				foreach my $sample (keys %{$sample_of_menthod1{$method}{$group}})
				{
					if    (!exists $sample_id_exp{$sample}{$diffgene_id}){print OUT"\t-4";}
					else  {$exp_out=log($sample_id_exp{$sample}{$diffgene_id})/log(10);print OUT"\t$exp_out";}
				}
			}
		print OUT "\n";
		}
		close IN;
		close OUT;
		open  my $fh_rcode,">$outdir/$prefix.$suffix.pheatmap.R" or die $!;
		print $fh_rcode <<CMD;
library("pheatmap")
data <- read.table("$outdir/$prefix.$suffix.prepare4pheatmap-plot",head=T,row.names=1)
aa <- head(data,1)
a <- t(aa)
sample <- row.names(a)[-1]
CMD
if (!defined $symbol){
print $fh_rcode <<CMD;
height <- 4
CMD
}
else {
print $fh_rcode <<CMD;
height <- length(row.names(data))/20
CMD
}
print $fh_rcode <<CMD;
if (length(data) >10){
width <- 2 + length(data)*0.2
}else{
width <- 4
}
group <- data.frame(Var1 = c($group_array_r)) 
rownames(group) = paste(sample)
colnames(group) = paste("Group")

Up_Down <- data.frame(data[,1]) 
rownames(Up_Down) = paste(row.names(data))
colnames(Up_Down) = paste("Up_Down")

################################################
###you can custom colors,or use default in R.###
################################################
ann_colors = list (Group = c("$group[0]"="#BEBADA","$group[1]"="#FDC086"),Up_Down=c("Up"="#E41A1C","Down"="#377EB8"))

pheatmap(
	data[,-1],
        annotation_col = group,
        annotation_row = Up_Down,
        annotation_names_col =F,
	annotation_names_row =F,
        annotation_colors = ann_colors,##group colors
        fontsize_col=9,
        fontsize_row=2.7,
		border_color = "transparent",
#        color = colorRampPalette(c("#FFFFD9","yellow2","orange","orangered2","gray10"))(200), #####colors and gradient for pheatmap 
	color = colorRampPalette(c("#1086BB","#A2D4E5","#F0F4D7","#FFD694","#FD432C"))(200), #####colors and gradient for pheatmap
main="Pheatmap for $vs",
CMD
if (!defined $symbol){
print $fh_rcode <<CMD;
	show_rownames = F,
CMD
}
print $fh_rcode <<CMD;
	file = '$outdir/$prefix.$suffix.pheatmap-plot.pdf',
	width=width,
	height=height
)
dev.off()
CMD
		#system("export R_LIBS=$Bin/3.2:$R_LIBS && $Rscript $outdir/$prefix.$suffix.pheatmap.R 2>/dev/null");# && rm $outdir/$prefix.$suffix.pheatmap-plot.R; ###licong
		system("export R_LIBS=$Bin/3.2:$R_LIBS && $Rscript $outdir/$prefix.$suffix.pheatmap.R");#chenweitian
                system("$Convert -density 300 -resize 30% $outdir/$prefix.$suffix.pheatmap-plot.pdf $outdir/$prefix.$suffix.pheatmap-plot.png");
		}
	}
}
print "done\n";
