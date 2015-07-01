#!/usr/bin/perl
# 
# Author: Long H. Do (long.h.do@gmail.com)
# 
# This script takes as input .diff files produced by CuffDiff to produce plots that search for GEDDs
# Requires R and the Bioconductor preprocessCore library
#
# usage: ./thisscript cuffdiff1 cuffdiff2 cuffdiff3
#
###########################################################################################
use strict;
use POSIX;

#takes cuffdiff .diff files and plots lowess smoothen lines across the chromosomes

my @y_title;# = ('Ts65Dn Fibroblast**','Ts65Dn Fibroblast*');
my @rho_order;# = ('1_3','2_3','1_2');
my @selected_chr = (); #shows all chromosomes
#my @selected_chr =(10,16); #plots only chromosomes 10 and 16

#########change this to your organism of choice to lable the chromosomes
my $box_prefix = 'Chr';
#my $box_prefix = 'HSA';

my $correlation_method = 'spearman';

my $pdf_filename = 'spearman.pdf';
my $pdf_width = 10;
my $pdf_height = 10;
my $pdf_column = 4;
my $pdf_row = 6;

my $plot_pointsize = .3;
my $plot_pointtype =  19;

###COLOR OPTIONS######
my %points_color = (1=>'lightpink',2=>'lightblue',3=>'lightgray');
my %lines_color = (1=>'red',2=>'blue',3=>'black');

my $box_title_size=.7;
my $box_title_pos = 0;

my $axis_label_font_size = .7;
my $axis_tick_size = .01;
my $axis_color = '#989898';

my $ylim = '-2,2';

my $rho_text_pos = '0,-1.7';
my $rho_font_size = 1.5;


my $plot_f = .1;
my $lowess_f = .3; #smoother span
my $lowess_iter = 1000; #iterations

my $tmpdir = '.tmp';

if(!-d $tmpdir){
	mkdir($tmpdir) || die "Error creating directory: $tmpdir\n"; 
}

open (RSCRIPT,"|/usr/bin/Rscript -") || die "$!\n";
#open (RSCRIPT,">rscript.txt");

print RSCRIPT '#! /usr/bin/Rscript --vanilla --default-packages=utils' . "\n";

_print_plots(_get_spearman(\@ARGV));

sub _print_plots{
	my ($spearman,$data) = @_;

	print RSCRIPT 'pdf("'.$pdf_filename.'",width='.$pdf_width.',height='.$pdf_height.');par(mfrow=c('.$pdf_column.','.$pdf_row.'),oma = c(4,7,0,0) + .1, mar = c(1,0,1,1) + .1)' . "\n";

	my @sorted;
	if($selected_chr[0]){
		@sorted = @selected_chr;
	}
	else{
		my @chr = keys(%{$spearman});
		@sorted = sort{$a <=> $b}@chr;
		if($sorted[0] eq "X"){
			shift(@sorted);
			push @sorted,"X";
		}
	}

	my (%uniq_num,%uniq_pair);

	for my $chr (@sorted){
		my @rho;
		
		my $n=1;
		for my $rho (@{$spearman->{$chr}}){
			my ($i,$j) = split '_', $rho;
			$uniq_pair{$rho}++;
			$uniq_num{$i}++;
			$uniq_num{$j}++;
			push @rho,'rho["'.$i.','.$j.'"] == .(round(rho_'."$chr\_$rho".',2))';
			$n++;
		}

		for my $i (@{$data->{$chr}}){
			if($i == 1){
				print RSCRIPT 'write("Plotting Chr' .$chr. '...",stderr()); log2data = log2(data.'."$chr".'[,'.($i*2).'] / data.'."$chr".'[,'.($i*2 - 1).']); plot(log2data,cex='.$plot_pointsize.',pch='.$plot_pointtype.',col="'.$points_color{$i}.'",ylim=c('.$ylim.'), xlab="",ylab="",axes=FALSE); box(lwd=.5,col="'.$axis_color.'");axis(side=1,lwd=.5,cex.axis='.$axis_label_font_size.', col="'.$axis_color.'",tck='.$axis_tick_size.',las=2,labels=NA);axis(side = 1, cex.axis='.$axis_label_font_size.', lwd = 0, line = -1.1, las = 1); axis(side=2,at=seq('.$ylim.',by=1),lwd=.5,cex.axis='.$axis_label_font_size.',col="'.$axis_color.'",tck='.$axis_tick_size.',las=3,labels=NA); axis(side = 2, at=seq('.$ylim.',by=1),cex.axis='.$axis_label_font_size.', lwd = 0, line = -.9, las = 2);lines(lowess(log2data,f='.$plot_f.',iter='.$lowess_iter.'),lwd=.4,col="'.$lines_color{$i}.'");' . "\n";
			}
			else{
				print RSCRIPT 'log2data = log2(data.'."$chr".'[,'.($i*2).'] / data.'."$chr".'[,'.($i*2 - 1).']) ; points(log2data,cex='.$plot_pointsize.',pch='.$plot_pointtype.',col="'.$points_color{$i}.'");lines(lowess(log2data,f='.$plot_f.',iter='.$lowess_iter.'),lwd=.4,col="'.$lines_color{$i}.'");' . "\n";
			}
		}
		print RSCRIPT 'mtext(side = 3, "'. $box_prefix . $chr.'", cex='.$box_title_size.',line ='.$box_title_pos.');text('.$rho_text_pos.',adj=0,cex='.$rho_font_size.',labels=bquote(' . join('*","~',@rho) . '));' . "\n";
	}

	print RSCRIPT 'mtext("Genes (chromosomal order)",side=1,cex=1,outer=TRUE,padj=1);';
	my $i = 1;
	my $padj = -.85;
	for my $y_title(@y_title){
		print RSCRIPT 'mtext(expression(paste(log[2],"[FC] ' . $y_title . '")),side=2,las=0,cex=1,adj=0,padj='. $padj .',outer=TRUE,at=0.43,col="'.$lines_color{$i}.'");';
		$padj -= 1.65;
		$i++;
	}

	unless( $selected_chr[0]){
		for my $num (keys %uniq_num){
			my @lowess_all;
			for my $chr (@sorted){
				push @lowess_all, 'lowess.'.$chr.'_'.$num . '$y';
			}
			print RSCRIPT "lowess.all_$num = c(" . join(",",@lowess_all) . ');';
		}

		for my $pair(sort keys %uniq_pair){
			my ($i,$j) = split '_',$pair;
			print RSCRIPT 'rho.all_' . "$i\_$j". '= cor(lowess.' . "all\_$i" . ',lowess.' . "all\_$j" . ',method="'.$correlation_method.'");' . "\n";
			print RSCRIPT "print('$i\_$j');print(rho.all_" ."$i\_$j);\n";
		}
	}

	close RSCRIPT;
	if(-d $tmpdir){
		system("rm -rf $tmpdir");
	}
}

sub _makeFC_chr{
	my ($files) = @_;
	
	my %selected;
	if($selected_chr[0]){
		for my $chr(@selected_chr){
			$selected{$chr}++;
		}
	}

	my %fpkm_chr;
	my @header;
	my $f=1;
	foreach my $file (@{$files}){
		open (FILE,"$file") || die "$!\n";
		my $n=1;
		my $flag=1;
		while (<FILE>){
			chomp;
			my @col = split;
			next unless ($col[7] >= .1); #only keep if fpkm > .1
			my($chr,$position) = split ':',$col[3];
			next unless ($chr && $position);
			$chr =~ s/chr//;
			next unless $chr =~ /^\d+$|^X$/;
			if($selected_chr[0]){
				next unless($selected{$chr});
			}

			my ($start,$stop) = split '-',$position;
			my $name = $col[4] . '__'. $col[5];
			if(!defined $fpkm_chr{$chr}{$start}{$f}){
				if($stop - $start > 0){
					push @{$fpkm_chr{$chr}{$start}{$f}},\@col;
				}
				else{
					push @{$fpkm_chr{$chr}{$stop}{$f}},\@col;
				}
			}

			if($n == 1){
				push @header,$name;
				push @y_title,$col[5] . " \($f\)";
			}
			$n++;
        }
    	close FILE;
		$f++;
	}

	##modify row/columns based on total chromosomes
	my $chr_count = scalar(keys(%fpkm_chr));
	if($chr_count <=6){
		$pdf_row = 2;
		$plot_pointsize = .2;
	}
	else{
		$pdf_row = 4;
		if($f>3){
			$rho_font_size = .9;
		}
	}
	$pdf_column = ceil($chr_count/$pdf_row);

	return(\%fpkm_chr,\@header);
}

sub _get_spearman{
	my ($files) = @_;
	my (%spearman,%data);

	my $pairs = scalar(@{$files});
	print RSCRIPT 'library(preprocessCore)' . "\n";

	my ($fpkm_chr,$header) = _makeFC_chr($files);

	my $header_string;
	for my $exp(@{$header}){
		my @e = split '__',$exp;
		$header_string .= "\t" . join("\t",@e);
	}
	for my $chr(keys %{$fpkm_chr}){
		next unless $chr =~ /^\d+$|^X$/;
		my $i = 1;
    	open (CHRFPKM,">$tmpdir/$chr\_fpkm.txt") || die "$!: Unable to create chr fpkm file\n";
		print CHRFPKM "gene\tlocus" . $header_string . "\n";
		my @position = keys(%{$fpkm_chr->{$chr}});
		my @sorted = sort{$a <=> $b}@position;
		for my $pos (@sorted){
			if( scalar(keys %{$fpkm_chr->{$chr}{$pos}}) == $pairs){
				my $n=1;
				for my $name (@{$header}){
					my @col = @{$fpkm_chr->{$chr}{$pos}{$n}};
					if($n==1){
						print CHRFPKM $col[0]->[2] . "\t" . $col[0]->[3];
					}
					for my $col(@col){
						print CHRFPKM "\t" . $col->[7] . "\t" . $col->[8];
					}
					$n++;
				}
				print CHRFPKM "\n";
			}
		}
		close CHRFPKM;
		print RSCRIPT 'write("Calculating lowess for Chr' .$chr. '...",stderr());' . 'chr_'.$chr.'_fpkm <- read.table("'. "$tmpdir\/$chr\_fpkm.txt" .'",header=T);' . "\n";
		#print 'write("Calculating lowess for Chr' .$chr. '...",stderr());' . 'chr_'.$chr.'_fpkm <- read.table("'. "$tmpdir\/$chr\_fpkm.txt" .'",header=T);' . "\n";

		print RSCRIPT 'data.'."$chr".' = normalize.quantiles(as.matrix('. "chr_$chr\_fpkm" .'[,3:'. (scalar(@{$files}) * 2 + 2) .']));' . "\n";
		#print 'data.'."$chr".' = normalize.quantiles(as.matrix('. "chr_$chr\_fpkm" .'[,3:'. (scalar(@{$files}) * 2 + 2) .']));' . "\n";

		for (my $i=1;$i<=$pairs;$i++){
			print RSCRIPT 'lowess.'."$chr\_$i".' = lowess(log2(data.'."$chr".'[,'. ($i*2) .']/data.'."$chr".'[,'. ($i*2 - 1).']),f='.$lowess_f.',iter='.$lowess_iter.');' . "\n";
			#print 'lowess.'."$chr\_$i".' = lowess(log2(data.'."$chr".'[,'. ($i*2) .']/data.'."$chr".'[,'. ($i*2 - 1).']),f='.$lowess_f.',iter='.$lowess_iter.');' . "\n";
			push @{$data{$chr}},$i;
		}

		my %combo;
		for (my $i=1;$i<=$pairs;$i++){
			for (my $j=1;$j<=$pairs;$j++){
				if($i<$j){
					$combo{"$i\_$j"}++;
				}
			}
		}

		my @order = sort keys %combo;

		my $n =1;
		for my $order(@order){
			my ($i,$j) = split '_',$order;
			print RSCRIPT 'rho_' ."$chr\_$i\_$j". '= cor(lowess.' . "$chr\_$i" . '$y,lowess.' . "$chr\_$j" . '$y,method="'.$correlation_method.'");' . "\n";
			#print 'rho_' ."$chr\_$i\_$j". '= cor(lowess.' . "$chr\_$i" . '$y,lowess.' . "$chr\_$j" . '$y,method="'.$correlation_method.'");' . "\n";
			push @{$spearman{$chr}},"$i\_$j";
			$n++;
		}
	}
	return (\%spearman,\%data);
}
