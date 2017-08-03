use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin '$Bin';
use Cwd 'abs_path';
use File::Basename;
my ($conf, $outdir, $help,$Rscript,$Convert,$zero);
GetOptions(
	"conf|c:s" => \$conf,
	"out|o:s"  => \$outdir,
	"Rpath:s" => \$Rscript,
	"convert:s" => \$Convert,
	"display_zero" => \$zero,
	"help|h:s" => \$help
);
if(!$conf || !$outdir || $help)
{
	die <<USAGE;
==================================================================================
Options:
* -conf    configure file of DiffVenn,consist of group_name,method_name,output file name,and input dir(diff result dir)
* -outdir  result directory
  -Rpath   path of Rscript
  -convert path of convert
  -display_zero  display "0" for zero values if 0, display nothing if undef
  -help    help information
E.g.:
perl $0 -conf /ifshk7/BC_RES/RD/USER/licong/RNA_RNAref/RNA_RNAref_2016a/GeneDiffVenn/GeneDiffVenn.conf  -out /ifshk7/BC_RES/RD/USER/licong/RNA_RNAref/RNA_RNAref_2016a/GeneDiffVenn/outdir
USAGE
}

$Rscript    ||= "$Bin/../software/Rscript";
$Convert    ||= "$Bin/../software/convert";
$zero       ||= 0;

my %para;
open CONF,$conf or die $!;
while (<CONF>){
next if (/^\s*$/ || /^\s*#/ || !/=/);
chomp;s/^([^#]+)#.*$/$1/;
my ($key,$value) = split /=/,$_,2;
$key =~ s/^\s*(.*?)\s*$/$1/;
$value =~ s/^\s*(.*?)\s*$/$1/;
next if ($value =~ /^\s*$/);
#print "$key\t$value\n";
$para{$key} = $value;
}


my (@group_name,$group_name,@plot_group,$plot_group,@every_res,$res,@resu,$venn_cnt,%hash,$file_name);
my $dault_name = 1;
if (defined $para{'Venn_VS_of_Method'})
{
	@plot_group = split /\s+/,$para{'Venn_VS_of_Method'};
	foreach $plot_group (@plot_group)
	{
		$dault_name = $dault_name + 1;
		if ($plot_group =~/\\\\|\/\//){
			@group_name = split /\\\\|\/\//,$plot_group;
			$file_name = $group_name[0];
			@every_res = split /,/,$group_name[1];
		}
		else {
			@every_res = split /,/,$plot_group;
			$file_name = $dault_name;
		}
		$venn_cnt = 0;
		foreach $res (@every_res){$venn_cnt ++;}
		foreach $res (@every_res)
       		{
				@resu = split /:/,$res;
				$hash{$file_name}{$venn_cnt}{$res}{$resu[1]}{$resu[0]} = 1;
       		}
	}
}


else {die "error in config:missing Venn_VS_of_Method";}
#print Dumper \%hash;
my(%hash_regulation,%hash_res,$venn_cnt,$res,$method,$vs);
if (defined %hash)
{
	my ($group_name,$venn_cnt,@diffline);
	foreach $group_name (keys %hash)
	{
		foreach $venn_cnt (keys %{$hash{$group_name}})
		{
			foreach $res (keys %{$hash{$group_name}{$venn_cnt}})
			{
				foreach $method (keys %{$hash{$group_name}{$venn_cnt}{$res}})
				{
					foreach $vs (keys %{$hash{$group_name}{$venn_cnt}{$res}{$method}})
					{
						$vs =~s/&/-VS-/g;
						#print "$vs_a[0]\t$vs_a[1]\n";
						my $diff_file = $para{'Diff_Dir'}."/".$vs.".".$method."_Method.GeneDiffExpFilter.xls";
						open IN,$diff_file, or die $diff_file;
						<IN>;
						while(my $line=<IN>)
						{
							chomp $line;
							@diffline = split /\t/,$line;
							$hash_regulation{$group_name}{$venn_cnt}{$diffline[0]}{$res}=$diffline[6];
						}
					}
				}
				push @{$hash_res{$group_name}{$venn_cnt}},$res;
			}
		}
	}
}


if (defined %hash_regulation)
{
	foreach my $group_name (keys %hash_regulation)
	{
		print "$group_name\n";
		foreach my $venn_cnt (keys %{$hash_regulation{$group_name}})
		{
			if ($venn_cnt == 2)
			{
				my($m1,$m2,@m1_v,@m2_v,$m1_r,$m2_r);
				$m1 = @{$hash_res{$group_name}{$venn_cnt}}[0];
				$m2 = @{$hash_res{$group_name}{$venn_cnt}}[1];
				@m1_v = split /:/,$m1;
				@m2_v = split /:/,$m2;
				initialize($zero,my ($up_1,$up_2,$up_1_2));
				initialize($zero,my ($down_1,$down_2,$down_1_2));
				initialize(0,my ($m1_c,$m2_c));
				foreach my $gi (keys %{$hash_regulation{$group_name}{$venn_cnt}})
				{
					if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m1}){$m1_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m1};}else{$m1_r = undef;}
					if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m2}){$m2_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m2};}else{$m2_r = undef;}
					if (defined $m1_r){$m1_c ++;}
					if (defined $m2_r){$m2_c ++;}
					if (defined $m1_r and !defined $m2_r)
					{
						if ($m1_r =~/up/i ) {$up_1 ++;}
						elsif($m1_r =~/down/i) {$down_1 ++;}
					}
					if (!defined $m1_r and defined $m2_r)
					{
						if ($m2_r =~/up/i ) {$up_2 ++;}
						elsif($m2_r =~/down/i) {$down_2 ++;}
					}
					if (defined $m1_r and defined $m2_r)
					{
						if ($m1_r =~/up/i ) {$up_1_2 ++;}
						elsif($m1_r =~/down/i) {$down_1_2 ++;}
					}
				}	
				##$group_name =~s/&/-VS-/g;
				open  my $fh_rcode,">$outdir/Venn_$group_name.R" or die "$outdir/Venn_$group_name.R";
				print $fh_rcode <<CMD;
library(plotrix)
pdf(file="$outdir/Venn_$group_name.pdf",width=10,height=8)
color <- c("#8DA0CB","lightcoral")
color_transparent <- adjustcolor(color, alpha.f = 0.15) 
color_transparent1 <- adjustcolor(color, alpha.f = 1)
########################################
p_x   <- c(-13,0,13)
p_y   <- c(0,0,0)
p_lab <- c("A","B","AB")
up_x   <- p_x
up_y   <- p_y +1
up_lab <- c("$up_1","$up_1_2","$up_2")
down_x <- p_x
down_y <- p_y -1
down_lab <- c("$down_1","$down_1_2","$down_2")

title_x <- c(-16,16)
title_y <- c(-16,-16)
title_lab <- c("$m1_v[1]\\n$m1_v[0]\\n($m1_c)","$m2_v[1]\\n$m2_v[0]\\n($m2_c)")
########################################
par(mar=c(6,8,6,13.5)+0.1,xpd=TRUE)
plot(c(-18,18), c(-18,18), type="n",xaxt = "n", xlab="",ylab="",yaxt = "n", axes=F,main="")
draw.ellipse(c(-5,5), c(0,0), c(14,14), c(14,14),border=color_transparent1,
angle=c(45,45), lty=1,col = color_transparent,lwd = 2)
text (up_x,up_y,up_lab,cex=0.8,col="#E41A1C")
text (down_x,down_y,down_lab,cex=0.8,col="#377EB8")
text (title_x,title_y,title_lab,cex=1.2,col=color_transparent1)
legend(26,4,legend=paste("Up"),col="#E41A1C",pch=20,bty = "n", pt.cex=2.8, cex=1.5)
legend(26,0, legend=paste("Down"),col="#377EB8",pch=20,bty = "n", pt.cex=2.8, cex=1.5) 
dev.off()
CMD
			
		}
		elsif ($venn_cnt == 3)
		{
			my ($m1,$m2,$m3,$m1_r,$m2_r,$m3_r,@m1_v,@m2_v,@m3_v);
			$m1 = @{$hash_res{$group_name}{$venn_cnt}}[0];
			$m2 = @{$hash_res{$group_name}{$venn_cnt}}[1];
			$m3 = @{$hash_res{$group_name}{$venn_cnt}}[2];
			@m1_v = split /:/,$m1;
			@m2_v = split /:/,$m2;
			@m3_v = split /:/,$m3;
			initialize($zero,my ($up_1,$up_2,$up_3,$up_12,$up_13,$up_23,$up_123));
			initialize($zero,my ($down_1,$down_2,$down_3,$down_12,$down_13,$down_23,$down_123));
			initialize(0,my ($m1_c,$m2_c,$m3_c));
			foreach my $gi (keys %{$hash_regulation{$group_name}{$venn_cnt}})
			{	#print "$gi\n";
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m1}){$m1_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m1};}else{$m1_r = undef;}
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m2}){$m2_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m2};}else{$m2_r = undef;}
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m3}){$m3_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m3};}else{$m3_r = undef;}
				if (defined $m1_r){$m1_c ++;}
				if (defined $m2_r){$m2_c ++;}
				if (defined $m3_r){$m3_c ++;}
				#print "$m1_v\t$m2_v\t$m3_v\n";
				if (defined $m1_r and !defined $m2_r and !defined $m3_r){
					if ($m1_r =~/up/i ) {$up_1 ++;}  
					elsif($m1_r =~/down/i) {$down_1 ++;}
				}
				if (!defined $m1_r and defined $m2_r and !defined $m3_r){
					if ($m2_r =~/up/i ) {$up_2 ++;}
					elsif($m2_r =~/down/i) {$down_2 ++;}
				}
				if (!defined $m1_r and !defined $m2_r and defined $m3_r){
					if ($m3_r =~/up/i ) {$up_3 ++;}
					elsif($m3_r =~/down/i) {$down_3 ++;}
				}
				if (defined $m1_r and defined $m2_r and !defined $m3_r){
					if ($m1_r =~ /up/i ) {$up_12 ++;}
					elsif($m1_r =~ /down/i) {$down_12 ++;}
				}
				if (!defined $m1_r and defined $m2_r and defined $m3_r){	
					if ($m2_r =~ /up/i ) {$up_23 ++;}
					elsif($m2_r =~ /down/i) {$down_23 ++;}
				}
				if (defined $m1_r and !defined $m2_r and defined $m3_r){	
					if ($m3_r =~ /up/i ) {$up_13 ++;}
					elsif($m3_r =~ /down/i) {$down_13 ++;}
				}
				if (defined $m1_r and defined $m2_r and defined $m3_r){	
					if ($m1_r =~ /up/i ) {$up_123 ++;}
					elsif($m1_r =~ /down/i) {$down_123 ++;}
				}
			}
				##$group_name =~s/&/-VS-/g;
				open  my $fh_rcode,">$outdir/Venn_$group_name.R" or die $!;
				print $fh_rcode <<CMD;
library(plotrix)
pdf(file="$outdir/Venn_$group_name.pdf",width=10,height=8)
color <- c("#E41A1C","#377EB8","#FDB462")
color_transparent <- adjustcolor(color, alpha.f = 0.15) 
color_transparent1 <- adjustcolor(color, alpha.f = 1)
########################################
p_x   <- c(0,13,-13,  10.4,-10.4,0, 0)
p_y   <- c(13,-9,-9,  4,4,-13.5, -1)
p_lab <- c("A","B","C","AB","AC","BC","ABC")

up_x   <- p_x
up_y   <- p_y +0.8
up_lab <- c ("$up_1","$up_2","$up_3","$up_12","$up_13","$up_23","$up_123")
down_x <- p_x
down_y <- p_y -0.8
down_lab <- c ("$down_1","$down_2","$down_3","$down_12","$up_13","$down_23","$down_123")

title_x <- c(0,17,-17)
title_y <- c(19.8,-18,-17)
title_lab <- c("$m1_v[1]\\n$m1_v[0]\\n($m1_c)","$m2_v[1]\\n$m2_v[0]\\n($m2_c)","$m3_v[1]\\n$m3_v[0]\\n($m3_c)")
########################################
par(mar=c(6,8,6,13.5)+0.1,xpd=TRUE)
plot(c(-18,18), c(-18,18), type="n",,xaxt = "n", xlab="",ylab="",yaxt = "n", axes=F,main="")
draw.ellipse(c(0,4,-4), c(3.02,-3.912,-3.912), c(14,14,14), c(14,14,14),border=color_transparent1,
angle=c(60,120,0), lty=1,col = color_transparent,lwd = 2)
text (up_x,up_y,up_lab,cex=0.8,col="#E41A1C")
text (down_x,down_y,down_lab,cex=0.8,col="#377EB8")
text (title_x,title_y,title_lab,cex=1.2,col=color_transparent1)
legend(26,4,legend=paste("Up"),col="#E41A1C",pch=20,bty = "n", pt.cex=2.8, cex=1.5)
legend(26,0, legend=paste("Down"),col="#377EB8",pch=20,bty = "n", pt.cex=2.8, cex=1.5)
dev.off()
CMD
	
		}
		elsif ($venn_cnt == 4)
		{
			my ($m1,$m2,$m3,$m4,$m1_r,$m2_r,$m3_r,$m4_r,@m1_v,@m2_v,@m3_v,@m4_v);
			$m1 = @{$hash_res{$group_name}{$venn_cnt}}[0];
           	$m2 = @{$hash_res{$group_name}{$venn_cnt}}[1];
            $m3 = @{$hash_res{$group_name}{$venn_cnt}}[2];
			$m4 = @{$hash_res{$group_name}{$venn_cnt}}[3];
			@m1_v = split /:/,$m1;
			@m2_v = split /:/,$m2;
			@m3_v = split /:/,$m3;
			@m4_v = split /:/,$m4;
			initialize($zero,my ($up_1,$up_2,$up_3,$up_4,  $up_12,$up_13,$up_14,$up_23,$up_24,$up_34, $up_123,$up_124,$up_134,$up_234, $up_1234));
			initialize($zero,my ($down_1,$down_2,$down_3,$down_4,  $down_12,$down_13,$down_14,$down_23,$down_24,$down_34, $down_123,$down_124,$down_134,$down_234, $down_1234));
			initialize(0,my ($m1_c,$m2_c,$m3_c,$m4_c));
			foreach my $gi (keys %{$hash_regulation{$group_name}{$venn_cnt}})
			{
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m1}){$m1_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m1};}else{$m1_r = undef;}
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m2}){$m2_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m2};}else{$m2_r = undef;}
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m3}){$m3_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m3};}else{$m3_r = undef;}
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m4}){$m4_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m4};}else{$m4_r = undef;}
				if (defined $m1_r){$m1_c ++;}
				if (defined $m2_r){$m2_c ++;}
				if (defined $m3_r){$m3_c ++;}
				if (defined $m4_r){$m4_c ++;}
				if (defined $m1_r and !defined $m2_r and !defined $m3_r and !defined $m4_r){
					if ($m1_r =~/up/i ) {$up_1 ++;}  
					elsif($m1_r =~/down/i) {$down_1 ++;}
				}
				if (!defined $m1_r and defined $m2_r and !defined $m3_r and !defined $m4_r){
					if ($m2_r =~/up/i ) {$up_2 ++;}           
                                        elsif($m2_r =~/down/i) {$down_2 ++;}
				}
				if (!defined $m1_r and !defined $m2_r and defined $m3_r and !defined $m4_r){
					if ($m3_r =~/up/i ) {$up_3 ++;}
					elsif($m3_r =~/down/i) {$down_3 ++;}
				}
				if (!defined $m1_r and !defined $m2_r and !defined $m3_r and defined $m4_r){
					if ($m4_r =~/up/i ) {$up_4 ++;}
					elsif($m4_r =~/down/i) {$down_4 ++;}
				}
				if (defined $m1_r and defined $m2_r and !defined $m3_r and !defined $m4_r){
					if ($m1_r =~/up/i ) {$up_12 ++;}
					elsif($m1_r =~/down/i) {$down_12 ++;}
				}
				if (defined $m1_r and !defined $m2_r and defined $m3_r and !defined $m4_r){
					if ($m2_r =~/up/i ) {$up_13 ++;}
					elsif($m2_r =~/down/i) {$down_13 ++;}
				}
				if (defined $m1_r and !defined $m2_r and !defined $m3_r and defined $m4_r){
					if ($m2_r =~/up/i ) {$up_14 ++;}
					elsif($m2_r =~/down/i) {$down_14 ++;}
				}
				if (!defined $m1_r and defined $m2_r and defined $m3_r and !defined $m4_r){
					if ($m3_r =~/up/i ) {$up_23 ++;}
					elsif($m3_r =~/down/i) {$down_23 ++;}
				}
				if (!defined $m1_r and defined $m2_r and !defined $m3_r and defined $m4_r){
					if ($m4_r =~/up/i ) {$up_24 ++;}
					elsif($m4_r =~/down/i) {$down_24 ++;}
				}
				if (!defined $m1_r and !defined $m2_r and defined $m3_r and defined $m4_r){
					if ($m1_r =~/up/i ) {$up_34 ++;}
					elsif($m1_r =~/down/i) {$down_34 ++;}
				}
                if (defined $m1_r and defined $m2_r and defined $m3_r and !defined $m4_r){
					if ($m2_r =~/up/i ) {$up_123 ++;}
					elsif($m2_r =~/down/i) {$down_123 ++;}
				}
				if (defined $m1_r and defined $m2_r and !defined $m3_r and defined $m4_r){
					if ($m2_r =~/up/i ) {$up_124 ++;}
					elsif($m2_r =~/down/i) {$down_124 ++;}
				}
				if (defined $m1_r and !defined $m2_r and defined $m3_r and defined $m4_r){
					if ($m2_r =~/up/i ) {$up_134 ++;}
					elsif($m2_r =~/down/i) {$down_134 ++;}
				}
                if (!defined $m1_r and defined $m2_r and defined $m3_r and defined $m4_r){
					if ($m3_r =~/up/i ) {$up_234 ++;}
					elsif($m3_r =~/down/i) {$down_234 ++;}
				}
				if (defined $m1_r and defined $m2_r and defined $m3_r and defined $m4_r){
					if ($m1_r =~/up/i ) {$up_1234 ++;}
					elsif($m1_r =~/down/i) {$down_1234 ++;}
				}
			}
			open  my $fh_rcode,">$outdir/Venn_$group_name.R" or die $!;
			print $fh_rcode <<CMD;
library(plotrix)
pdf(file="$outdir/Venn_$group_name.pdf",width=10,height=8)
color <- c("#E41A1C","#377EB8","#FDB462","#4DAF4A") 
color_transparent <- adjustcolor(color, alpha.f = 0.15) 
color_transparent1 <- adjustcolor(color, alpha.f = 1)

p_x   <- c(-15,15,-7.3,7.3,  0,-10,-9,10,9,0,  4.9,-4.9,-6.5,6.5 ,0)
p_y   <- c(2,2,11,11,  -12,6.3,-4,6.3,-4,8.8,  -7.4,-7.4,2.6,2.6 ,-2)
p_lab <- c ("A","B","C","D","AB","AC","AD","BC","BD","DE",  "ABC","ABD","ACD","BCD",  "ABCD")

up_x   <- p_x
up_y   <- p_y +0.8
up_lab <- c("$up_1","$up_2","$up_3","$up_4","$up_12","$up_13","$up_14","$up_23","$up_24","$up_34","$up_123","$up_124","$up_134","$up_234","$up_1234")
down_x <- p_x
down_y <- p_y -0.8
down_lab <- c("$down_1","$down_2","$down_3","$down_4","$down_12","$down_13","$down_14","$down_23","$down_24","$down_34","$down_123","$down_124","$down_134","$down_234","$down_1234")

title_x <- c(-18,18,-15,15)
title_y <- c(-16,-16,16,16)
title_lab <- c("$m1_v[1]\\n$m1_v[0]\\n($m1_c)","$m2_v[1]\\n$m2_v[0]\\n($m2_c)","$m3_v[1]\\n$m3_v[0]\\n($m3_c)","$m4_v[1]\\n$m4_v[0]\\n($m4_c)")
#########################################
par(mar=c(6,8,6,13.5)+0.1,xpd=TRUE)
plot(c(-18,18), c(-18,18), type="n",xaxt = "n", xlab="",ylab="",yaxt = "n", axes=F,main="")
draw.ellipse(c(-6,6,0,0), c(-4,-4,2,2), c(16,16,14,14), c(9.89,9.89,8.65,8.65),border=color_transparent1,
angle=c(-40,40,-45,45), lty=1,col = color_transparent,lwd = 2)
text (up_x,up_y,up_lab,cex=0.8,col="#E41A1C")
text (down_x,down_y,down_lab,cex=0.8,col="#377EB8")
text (title_x,title_y,title_lab,cex=1.2,col=color_transparent1)
legend(26,4,legend=paste("Up"),col="#E41A1C",pch=20,bty = "n", pt.cex=2.8, cex=1.5)
legend(26,0, legend=paste("Down"),col="#377EB8",pch=20,bty = "n", pt.cex=2.8, cex=1.5) 
dev.off()
CMD
		}
		elsif ($venn_cnt == 5)
		{
			my ($m1,$m2,$m3,$m4,$m5,$m1_r,$m2_r,$m3_r,$m4_r,$m5_r,@m1_v,@m2_v,@m3_v,@m4_v,@m5_v);	
			$m1 = @{$hash_res{$group_name}{$venn_cnt}}[0];
			$m2 = @{$hash_res{$group_name}{$venn_cnt}}[1];
			$m3 = @{$hash_res{$group_name}{$venn_cnt}}[2];
			$m4 = @{$hash_res{$group_name}{$venn_cnt}}[3];
			$m5 = @{$hash_res{$group_name}{$venn_cnt}}[4];
			@m1_v = split /:/,$m1;
			@m2_v = split /:/,$m2;
			@m3_v = split /:/,$m3;
			@m4_v = split /:/,$m4;
			@m5_v = split /:/,$m5;
			#initialize(0,my ($up_1,$up_2,$up_3,$up_4,$up_5,$up_12,$up_13,$up_14,$up_15,$up_23,$up_24,$up_25,$up_34,$up_35,$up_45,$up_123,$up_124,$up_125,$up_134,$up_135,$up_145,$up_234,$up_235,$up_245,$up_345,$up_1234,$up_1235,$up_1245,$up_1345,$up_2345,$up_12345));
			#initialize(0,my ($down_1,$down_2,$down_3,$down_4,$down_5,$down_12,$down_13,$down_14,$down_15,$down_23,$down_24,$down_25,$down_34,$down_35,$down_45,$down_123,$down_124,$down_125,$down_134,$down_135,$down_145,$down_234,$down_235,$down_245,$down_345,$down_1234,$down_1235,$down_1245,$down_1345,$down_2345,$down_12345));
			initialize(0,my ($m1_c,$m2_c,$m3_c,$m4_c,$m5_c));
			initialize($zero,my ($up_1,$up_2,$up_3,$up_4,$up_5,$up_12,$up_13,$up_14,$up_15,$up_23,$up_24,$up_25,$up_34,$up_35,$up_45,$up_123,$up_124,$up_125,$up_134,$up_135,$up_145,$up_234,$up_235,$up_245,$up_345,$up_1234,$up_1235,$up_1245,$up_1345,$up_2345,$up_12345));
			initialize($zero,my ($down_1,$down_2,$down_3,$down_4,$down_5,$down_12,$down_13,$down_14,$down_15,$down_23,$down_24,$down_25,$down_34,$down_35,$down_45,$down_123,$down_124,$down_125,$down_134,$down_135,$down_145,$down_234,$down_235,$down_245,$down_345,$down_1234,$down_1235,$down_1245,$down_1345,$down_2345,$down_12345));
			foreach my $gi (keys %{$hash_regulation{$group_name}{$venn_cnt}})
			{
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m1}){$m1_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m1};}else{$m1_r = undef;}
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m2}){$m2_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m2};}else{$m2_r = undef;}
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m3}){$m3_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m3};}else{$m3_r = undef;}
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m4}){$m4_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m4};}else{$m4_r = undef;}
				if (exists $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m5}){$m5_r = $hash_regulation{$group_name}{$venn_cnt}{$gi}{$m5};}else{$m5_r = undef;}
				if (defined $m1_r){$m1_c ++;}
				if (defined $m2_r){$m2_c ++;}
				if (defined $m3_r){$m3_c ++;}
				if (defined $m4_r){$m4_c ++;}
				if (defined $m5_r){$m5_c ++;}
				if (defined $m1_r and !defined $m2_r and !defined $m3_r and !defined $m4_r and !defined $m5_r){
					if ($m1_r =~/up/i ) {$up_1 ++;}  
					elsif($m1_r =~/down/i) {$down_1 ++;}
				}
				if (!defined $m1_r and defined $m2_r and !defined $m3_r and !defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_2 ++;}
					elsif($m2_r =~/down/i) {$down_2 ++;}
				}
				if (!defined $m1_r and !defined $m2_r and defined $m3_r and !defined $m4_r and !defined $m5_r){
					if ($m3_r =~/up/i ) {$up_3 ++;}
					elsif($m3_r =~/down/i) {$down_3 ++;}
				}
				if (!defined $m1_r and !defined $m2_r and !defined $m3_r and defined $m4_r and !defined $m5_r){
					if ($m4_r =~/up/i ) {$up_4 ++;}
					elsif($m4_r =~/down/i) {$down_4 ++;}
				}
				if (!defined $m1_r and !defined $m2_r and !defined $m3_r and !defined $m4_r and defined $m5_r){
					if ($m5_r =~/up/i ) {$up_5 ++;}
					elsif($m5_r =~/down/i) {$down_5 ++;}
				}
				if (defined $m1_r and defined $m2_r and !defined $m3_r and !defined $m4_r and !defined $m5_r){
					if ($m1_r =~/up/i ) {$up_12 ++;}
					elsif($m1_r =~/down/i) {$down_12 ++;}
				}
				if (defined $m1_r and !defined $m2_r and defined $m3_r and !defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_13 ++;}
					elsif($m2_r =~/down/i) {$down_13 ++;}
				}
				if (defined $m1_r and !defined $m2_r and !defined $m3_r and defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_14 ++;}
					elsif($m2_r =~/down/i) {$down_14  ++;}
				}
				if (defined $m1_r and !defined $m2_r and !defined $m3_r and !defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_15 ++;}
					elsif($m2_r =~/down/i) {$down_15 ++;}
				}
				if (!defined $m1_r and defined $m2_r and defined $m3_r and !defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_23 ++;}
					elsif($m2_r =~/down/i) {$down_23  ++;}
				}
				if (!defined $m1_r and defined $m2_r and !defined $m3_r and defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_24 ++;}
					elsif($m2_r =~/down/i) {$down_24 ++;}
				}
				if (!defined $m1_r and defined $m2_r and !defined $m3_r and !defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_25 ++;}
					elsif($m2_r =~/down/i) {$down_25 ++;}
				}
				if (!defined $m1_r and !defined $m2_r and defined $m3_r and defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_34 ++;}
					elsif($m2_r =~/down/i) {$down_34 ++;}
				}
				if (!defined $m1_r and !defined $m2_r and defined $m3_r and !defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_35 ++;}
					elsif($m2_r =~/down/i) {$down_35 ++;}
				}
				if (!defined $m1_r and !defined $m2_r and !defined $m3_r and defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_45 ++;}
					elsif($m2_r =~/down/i) {$down_45 ++;}
				}
				if (defined $m1_r and defined $m2_r and defined $m3_r and !defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_123 ++;}
					elsif($m2_r =~/down/i) {$down_123 ++;}
				}
				if (defined $m1_r and defined $m2_r and !defined $m3_r and defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_124 ++;}
					elsif($m2_r =~/down/i) {$down_124 ++;}
				}
				if (defined $m1_r and defined $m2_r and !defined $m3_r and !defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_125 ++;}
					elsif($m2_r =~/down/i) {$down_125 ++;}
				}
				if (defined $m1_r and !defined $m2_r and defined $m3_r and defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_134 ++;}
					elsif($m2_r =~/down/i) {$down_134 ++;}
				}
				if (defined $m1_r and !defined $m2_r and defined $m3_r and !defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_135 ++;}
					elsif($m2_r =~/down/i) {$down_135 ++;}
				}
				if (defined $m1_r and !defined $m2_r and !defined $m3_r and defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_145 ++;}
					elsif($m2_r =~/down/i) {$down_145 ++;}
				}
				if (!defined $m1_r and defined $m2_r and defined $m3_r and defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_234 ++;}
					elsif($m2_r =~/down/i) {$down_234  ++;}
				}
				if (!defined $m1_r and defined $m2_r and defined $m3_r and !defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_235 ++;}
					elsif($m2_r =~/down/i) {$down_235 ++;}
				}
				if (!defined $m1_r and defined $m2_r and !defined $m3_r and defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_245 ++;}
					elsif($m2_r =~/down/i) {$down_345 ++;}
				}
				if (defined $m1_r and defined $m2_r and defined $m3_r and defined $m4_r and !defined $m5_r){
					if ($m2_r =~/up/i ) {$up_1234 ++;}
					elsif($m2_r =~/down/i) {$down_1234 ++;}
				}
				if (defined $m1_r and defined $m2_r and defined $m3_r and !defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_1235 ++;}
					elsif($m2_r =~/down/i) {$down_1235 ++;}
				}
				if (defined $m1_r and defined $m2_r and !defined $m3_r and defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_1245 ++;}
					elsif($m2_r =~/down/i) {$down_1245 ++;}
				}
				if (defined $m1_r and !defined $m2_r and defined $m3_r and defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_1345 ++;}
					elsif($m2_r =~/down/i) {$down_1345 ++;}
				}
				if (!defined $m1_r and defined $m2_r and defined $m3_r and defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_2345 ++;}
					elsif($m2_r =~/down/i) {$down_2345 ++;}
				}
				if (defined $m1_r and defined $m2_r and defined $m3_r and defined $m4_r and defined $m5_r){
					if ($m2_r =~/up/i ) {$up_12345 ++;}
					elsif($m2_r =~/down/i) {$down_12345 ++;}
				}
			}
			open  my $fh_rcode,">$outdir/Venn_$group_name.R" or die $!;
			print $fh_rcode <<CMD;
library(plotrix)
pdf(file="$outdir/Venn_$group_name.pdf",width=10,height=8)
color <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FDB462") 
color_transparent <- adjustcolor(color, alpha.f = 0.15) 
color_transparent1 <- adjustcolor(color, alpha.f = 1)
########################################################
p_x   <- c(-4.5,14,13,-6,-18,  7.2,7.2,0,-9,13.2,-7.5,12.05,1.5,-12,-12.5,  9,3.7,9.1,4,-10.5,-4,-3.8,12.2,-10,-11.5,  2.2,9.9,3.7,-7.8,-8.5,  0)
p_y   <- c(17,9.5,-10.5,-15.5,1,  11,-10.8,12.7,10,-3.4,-10.8,4,-13.5,3.5,-5,  -7.2,12.1,6.5,-11.7,7.4,10.3,-11,-0.1,-7.2,0.3,  -10.1,-1,9,6.3,-5.8,  0)
p_lab <- c ("A","B","C","D","E","AB","AC","AD","AE","BC","BD","BE","CD","CE","DE","ABC","ABD","ABE","ACD","ACE","ADE","BCD","BCE","BDE","CDE","ABCD","ABCE","ABDE","ACDE","BCDE","ABCDE")

up_x   <- p_x
up_y   <- p_y +0.6
up_lab <- c("$up_1","$up_2","$up_3","$up_4","$up_5","$up_12","$up_13","$up_14","$up_15","$up_23","$up_24","$up_25","$up_34","$up_35","$up_45","$up_123","$up_124","$up_125","$up_134","$up_134","$up_145","$up_234","$up_235","$up_245","$up_345","$up_1234","$up_1235","$up_1245","$up_1345","$up_2345","$up_12345")
down_x <- p_x
down_y <- p_y -0.6
down_lab <- c("$down_1","$down_2","$down_3","$down_4","$down_5","$down_12","$down_13","$down_14","$down_15","$down_23","$down_24","$down_25","$down_34","$down_35","$down_45","$down_123","$down_124","$down_125","$down_134","$down_134","$down_145","$down_234","$down_235","$down_245","$down_345","$down_1234","$down_1235","$down_1245","$down_1345","$down_2345","$down_12345")
############################################
title_x <-c(-16,17,22,-14,-22)
title_y <-c(21,18.5,-18,-23,-8)
title_lab <- c("$m1_v[1]\\n$m1_v[0]\\n($m1_c)","$m2_v[1]\\n$m2_v[0]\\n($m2_c)","$m3_v[1]\\n$m3_v[0]\\n($m3_c)","$m4_v[1]\\n$m4_v[0]\\n($m4_c)","$m5_v[1]\\n$m5_v[0]\\n($m5_c)")
########################################################
par(mar=c(6,8,6,13.5)+0.1,xpd=TRUE)
plot(c(-18,18), c(-18,18), type="n",xaxt = "n", xlab="",ylab="",yaxt = "n", axes=F,main="")
draw.ellipse(c(0,4,2.4,-2.4,-4), c(4.35,1.1,-3.55,-3.55,1.1), c(18,18,18,18,18), c(10.5,10.5,10.5,10.5,10.5),border=color_transparent1,
 angle=c(290,218,146,74,2), lty=1,col = color_transparent,lwd = 2)
text (up_x,up_y,up_lab,cex=0.8,col="#E41A1C")
text (down_x,down_y,down_lab,cex=0.8,col="#377EB8")
text (title_x,title_y,title_lab,cex=1.2,col=color_transparent1)
legend(26,4,legend=paste("Up"),col="#E41A1C",pch=20,bty = "n", pt.cex=2.8, cex=1.5)
legend(26,0, legend=paste("Down"),col="#377EB8",pch=20,bty = "n", pt.cex=2.8, cex=1.5)
dev.off() 
CMD
	}
}
`$Rscript $outdir/Venn_$group_name.R 2>/dev/null`;
system("$Convert -density 300 -resize 40% $outdir/Venn_$group_name.pdf $outdir/Venn_$group_name.png");
}
}


sub initialize
{
my $value = shift(@_);
  foreach my $elem (@_)
  {
      $elem = $value;
  }
}
