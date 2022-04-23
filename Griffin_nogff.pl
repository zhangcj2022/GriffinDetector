#!/usr/bin/perl -w
#Draw the figure GriffinDetector v1.1

my $ctl_file=$ARGV[0];
print "Please input the control file\n" unless $ctl_file;
my $skip_step2="N";
my $indenfity_doc="";

if (exists $ARGV[1]){
    if ($ARGV[1] eq "Y"){
	$skip_step2="Y";
	unless ($ARGV[2]){
	    print "Please input the indentified_chimerics path\n";
	    exit;
	}
	$indenfity_doc=$ARGV[2];
	unless (-d $indenfity_doc){
	    print "The indentified_chimerics path is not a directory\n";
	    exit;
	}
    }
}

my $flag_spe="F";
my $focus_species;
my @lines;
my %focus_species_name;
my %outgroup_species_name;
my $makeblastdb;
my $blastp;
my $skip;
my $outfile;
my @intree;
my $phylogeny;
my %protein_path;
my @intree_check;
my @group2;#这个参数是为了和原先的classfy_chimeric_genes.pl对应而使用的，里面将会包含$file4
###################################################################################################################
#==# reading control file

open (CTL,$ctl_file) or die "Cannot open the control file: $ctl_file, please check if the path is correct or not\n";
while (<CTL>){
    next if $_ =~/\#/;
    while ($_=~s/\s+//gi){}
    if ($_=~/SPECIES/ and $flag_spe eq "F"){
        $flag_spe="T";
        next;
    }
    chomp;
    if ($_=~/SPECIES/ and $flag_spe eq "T"){
        $flag_spe="F";
        next;
    }
    if ($flag_spe eq "T"){
        @lines=split/\:/,$_;
        $protein_path{$lines[0]}=$lines[1];
        push (@intree_check,$lines[0]);
        if (-d "$lines[1]"){
            print "Cannot find the fasta file of $lines[0], please check\n";
            exit;
        }
        if (-e "$lines[1]"){
            #print $lines[1],"\n";
        }else{
            print "Cannot find the fasta file of $lines[0], please check\n";
            exit;
        }

    }
    if ($_=~s/TREE://){
        @intree=&tree($_);
        $phylogeny=$_;
        while ($phylogeny=~s/ //g){}
        next;
    }
    if ($_=~/INGROUPS/){
        undef @lines;
        @lines=split/\:|\,/;
        $focus_species=@lines;
        for ($i=1;$i<$focus_species;$i++){
            $focus_species_name{$lines[$i]}++;
        }
    }
    if ($_=~/OUTGROUPS/){
        undef @lines;
        @lines=split/\:|\,/;
        $focus_species=@lines;
        for ($i=1;$i<$focus_species;$i++){
            $outgroup_species_name{$lines[$i]}++;
        }
    }
    if ($_=~s/Makeblastdb://){
        $makeblastdb=$_;
        if (-e "$makeblastdb"){}else{
            print "Cannot find the Makeblastdb program, please check\n";
            exit;
        }
    }
    if ($_=~s/Blastp://){
        $blastp=$_;
        if (-e "$blastp"){}else{
            print "Cannot find the BLASTP program, please check\n";
            exit;
        }
    }
    if ($_=~s/SKIP://){
        $skip=$_;
    }
    if ($_=~s/OUTFILE://){
        $outfile=$_;
    }
}
close CTL;

###################################################################################################################
#==# treat parameters
if ($skip eq "Y"){
    if ($outfile eq ""){
        print "Please input the OUTFILE option since the SKIP was set to Y; otherwise, please set the SKIP to N\n";
        exit;
    }
}

my @tmpa=sort @intree;
my @tmpb=sort @intree_check;
if (@tmpa ne @tmpb){
    print "The species number is not identical in SPECIES and TREE!\n";
    exit;
}
my %intree_hash;
for ($i=0;$i<=$#tmpa;$i++){
    $intree_hash{$tmpa[$i]}++;
    if ($tmpa[$i] ne $tmpb[$i]){
        print "The species name are not identical in SPECIES and TREE!\n";
        exit;
    }
}

foreach my $key (sort keys %focus_species_name){
    if (exists $outgroup_species_name{$key}){
        print "There're species in INGROUP is also presented in OUTGROUP\n";
        exit;
    }
    if (not exists $intree_hash{$key}){
        print "Species in INGROUP is not presented in TREE!\n";
        exit;
    }
}
foreach my $key (sort keys %outgroup_species_name){
    if (not exists $intree_hash{$key}){
        print "Species in OUTGROUP is not presented in TREE!\n";
        exit;
    }
}

if ($skip_step2 eq "Y"){
    goto Lable;
}
###################################################################################################################
#==# treat data
my $time_id=time();
use File::Copy 'move';
if ($skip eq "Y"){
    chdir $outfile;
    my @formated_files=("",".pal",".phd",".phi",".phr",".pin",".pog",".psd",".psi",".psq");
    my $tmp_file;
    foreach $key (sort keys %protein_path){
        my @tmp_dbfile=split/\//,$protein_path{$key};
        my $dbfile=$tmp_dbfile[$#tmp_dbfile];
        my $postfix;
        if (-d $key){
            chdir $key;
            foreach $postfix (@formated_files){
                $tmp_file=$dbfile.$postfix;
                unless (-e $tmp_file){
                    print "The data is insufficient (makeblastdb of species $key is not finished), please set the SKIP option into N and run again; or you can manually formated the database and run again.\n";
                    exit;
                }
            }
            chdir "..";
        }else{
            print "The data is insufficient (species $key was not formated), please set the SKIP option into N and run again;  or you can manually formated the database and run again.\n";
            exit;
        }
        foreach $tmp_file (sort keys %focus_species_name){
            $dbfile=$tmp_file."2".$key."_aa2aa_blast+.query";
            unless (-e $dbfile){
                print "The data is insufficient (blastp results between $tmp_file and species $key was not presented), please set the SKIP option into N and run again;  or you can manually run the blastp, make sure the outfile name should be $dbfile.\n";
                exit;
            }
        }
    }
}else{
    if (-d "blast_plus_db"){
        move "blast_plus_db", "blast_plus_db".$time_id or die ("back blast_plus_db failed: $!");
        print "Your previous DB was back up as blast_plus_db".$time_id,"\n";
        mkdir "blast_plus_db";
    }else{
        mkdir "blast_plus_db";
    }
    chdir "blast_plus_db";
    my $cmd;
    my $processor=1;
    my @child;
    my %db_hash;
    foreach $key (sort keys %protein_path){
        my @tmp_dbfile=split/\//,$protein_path{$key};
        my $dbfile=$tmp_dbfile[$#tmp_dbfile];
        my $tmp_path=$protein_path{$key};
        $tmp_path=~s/$dbfile//;
        $db_hash{$key}=$key."/".$dbfile;
    }
    foreach $key (sort keys %protein_path){
        mkdir "$key";
        chdir "$key";
        $cmd="cp $protein_path{$key} .";
        system $cmd;
        my @tmp_dbfile=split/\//,$protein_path{$key};
        my $dbfile=$tmp_dbfile[$#tmp_dbfile];

        $cmd="$makeblastdb -in $dbfile -hash_index -dbtype prot -parse_seqids -gi_mask";
        system $cmd;
        chdir "..";
    }

    $processor=1;
    foreach $key (sort keys %focus_species_name){
        foreach my $tmp_key (sort keys %protein_path){
            $cmd="$blastp -db $db_hash{$tmp_key} -query $db_hash{$key} -out $key"."2".$tmp_key."_aa2aa_blast+.query -evalue 0.00001 -outfmt 7 -max_target_seqs 10 -num_threads 10";
            print $cmd,"\n\n";
#           exit;
            if ($processor <8){
                $child[$processor]=fork();
                if ($child[$processor]==0){
                    exec $cmd;
                    exit 0;
                }
                $processor++;
            }else{
                for (my $i=1;$i<$processor;$i++){
                    waitpid($child[$i],0);
                }
                $processor=1;
                $child[$processor]=fork();
                if ($child[$processor]==0){
                    exec $cmd;
                    exit 0;
                }
                $processor++;
            }
        }
        for (my $i=1;$i<$processor;$i++){
            waitpid($child[$i],0);
        }
    }
    $outfile="blast_plus_db";
}
chdir "..";
if (-d "identified_chimerics"){
    move "identified_chimerics", "identified_chimerics".$time_id or die ("back identified_chimerics failed: $!");
    print "Your previous identified_chimerics was back up as identified_chimerics".$time_id,"\n";
    mkdir "identified_chimerics";
}else{
    mkdir "identified_chimerics";
}
chdir "identified_chimerics";
my %candidate_chimeric_gene_id;
my %chimeric_childs_ids_hash_array;
my %query_aa;
my %target_aa;
my %hash_scan;
my $check_file2;
foreach $key (sort keys %focus_species_name){
    mkdir "$key";
    chdir "$key";
    my $tmp_name;
    foreach $tmp_name (sort keys %protein_path){
        my $file1=$protein_path{$key};
        my $file2=$protein_path{$tmp_name};
        my $file3;
	$file3="../../".$outfile."/".$key."2".$tmp_name."_aa2aa_blast+.query";
        my $file4=$key."2".$tmp_name."candidate_chimerics_long_copy.ids";
        push (@group2,"$tmp_name:$key/$file4");
        my $file5=$key."2".$tmp_name."candidate_chimeric_short_copy.ids";
        $check_file2=$file2;
        undef %query_aa;
        undef %target_aa;
        undef %hash_scan;
        %query_aa=&read_fasta_file_new($file1,"space");
        %target_aa=&read_fasta_file_new($file2,"space");
        my $protein_hits=0;
        open (FILE,$file3) or die "cannot open file $file3\n";
        open (COPIES,">$file4") or die "cannot open file $file4\n";
        open (SCAN,">$file5") or die "cannot open file $file5\n";
        my $tmpid;
        my $hit_count;
        my %hit_count;
	my $zerohit;
        while (<FILE>){
            if ($_=~/BLASTP/){
                <FILE>;
                <FILE>;
                $zerohit=<FILE>;
                if ($zerohit =~/hits/){
		    #print "OK" if $zerohit eq "# 0 hits found\n";
                }else{
                    <FILE>;
                }
                if ($protein_hits ne "0"){
                    &scan_chimeric($tmpid);
                }
                $protein_hits=0;
                next;
            }elsif ($_=~/processed/){
		if ($zerohit eq "# 0 hits found\n"){
		}else{
		    &scan_chimeric($tmpid);
		}
		
	    }else{
                chomp;
                my @lines=split/\s+/,$_;
                $protein_hits++;
                $hit_count{$lines[0]}{$lines[1]}++;
                $hit_count=$hit_count{$lines[0]}{$lines[1]};
                if ($hit_count>1){
                    $hash_scan{$lines[0]}{$lines[1]}{q_start}=$lines[6] if $hash_scan{$lines[0]}{$lines[1]}{q_start}>$lines[6];
                    $hash_scan{$lines[0]}{$lines[1]}{q_stop}=$lines[7] if $hash_scan{$lines[0]}{$lines[1]}{q_stop}<$lines[7];
                    $hash_scan{$lines[0]}{$lines[1]}{t_start}=$lines[8] if $hash_scan{$lines[0]}{$lines[1]}{t_start}>$lines[8];
                    $hash_scan{$lines[0]}{$lines[1]}{t_stop}=$lines[9] if $hash_scan{$lines[0]}{$lines[1]}{t_stop}<$lines[9];
#                    $hash_scan{$lines[0]}{$lines[1]}{identity}.="+".$lines[2];目前的任何版本都没有使用这一数据
                }else{
                    $hash_scan{$lines[0]}{$lines[1]}{q_start}=$lines[6];
                    $hash_scan{$lines[0]}{$lines[1]}{q_stop}=$lines[7];
                    $hash_scan{$lines[0]}{$lines[1]}{t_start}=$lines[8];
                    $hash_scan{$lines[0]}{$lines[1]}{t_stop}=$lines[9];
                    $hash_scan{$lines[0]}{$lines[1]}{identity}=$lines[2];
                }
                $tmpid=$lines[0];

            }
        }

        close FILE;
        close COPIES;
        close SCAN;
    }
    chdir "..";
}

##############################################################################################################
#==# generate enough data so we can seperate
Lable:
    if ($skip_step2 eq "Y"){
	chdir $indenfity_doc;
}
my @group1;
foreach $key (sort keys %focus_species_name){
    my $tmp_name;
    my %Bfinal_chimeric_hash;
    foreach $tmp_name (sort keys %protein_path){
        my $file5=$key."/".$key."2".$tmp_name."candidate_chimeric_short_copy.ids";
	open (FILE,$file5) or die "cannot open file $file5\n";#是short copy ids file，$file4是long copy ids
        my %chimeric;
        my %info;
        while (<FILE>){
            chomp;
            @lines=split/\s+/,$_;
            $chimeric{$lines[0]}++;
            $info{$lines[0]}{$chimeric{$lines[0]}}=$_;
        }
        close FILE;
        my $tmpkey3;
	my $file6=$key."/".$key."2".$tmp_name."more_than_2_short_copy.ids";
	open (FILE,">$file6") or die "cannot open $file6\n";
        foreach $tmpkey3 (sort {$chimeric{$b}<=>$chimeric{$a}} keys %chimeric){#这个是在目前这个物种中所有的chimeric gene id
            if ($chimeric{$tmpkey3}>1){#把有两个或两个以上没有overlap的short copy取出来
		my @tmpinfo;
		my @newtmp_info;
		my %tmphash;
		my $ttmp;
		for ($i=1;$i<=$chimeric{$tmpkey3};$i++){

#问题的根源，当全部都没有“点”的时候，就是把short copy的id直接输入来检测了
#例如 
#chengjunzh1982@minion12$ head ../flybase_data/identified_chimerics/dmel/dmel2Aedescandidate_chimeric_short_copy.ids
#FBpp0289638     AAEL006790-PA   0.771   0.985   0.214   0.753   0.999
#就是直接把AAEL006790-PA输入了$tmphash,来看是否两个短copy来自于同一个transcript

#当长copy不存在点，而short copy存在点的时候，
#例如 
#chengjunzh1982@minion12$ head ../identified_chimerics_for_rapdb/jap/jap2barthiicandidate_chimeric_short_copy.ids
#Os12t0512000-01 Obart12g13520.1 0.451   1       0.549   0.004   1

#当两个copy都存在点的时候，那么事情就变成了
#chengjunzh1982@minion12$ head ../recalculate/identified_chimerics/jap/jap2barthiicandidate_chimeric_short_copy.ids
#LOC_Os01g01010.1        Obart02g29750.1 0.087   0.778   0.691   0.514   0.88
#
#假如长copy有点，而短copy没有点，那么，就是
#以下是杜撰的内容
#Obart12g13520.1        Os12t0512000-01 0.087   0.778   0.691   0.514   0.88

		    @tmpinfo=split/\.|\s+/,$info{$tmpkey3}{$i};
		    @newtmpinfo=split/\s+/,$info{$tmpkey3}{$i};		    
		    if ($newtmpinfo[0]=~/\./){
#这里,果蝇的版本将于水稻的不同
			$tmphash{$tmpinfo[2]}=$info{$tmpkey3}{$i} if not exists $tmphash{$tmpinfo[2]};

		    }else{
			$tmphash{$tmpinfo[1]}=$info{$tmpkey3}{$i} if not exists $tmphash{$tmpinfo[1]};
		    }
		}
		$ttmp=keys %tmphash;
		if ($ttmp>1){
		    my $tmpkey4;

		    foreach $tmpkey4 (sort keys %tmphash){
			print FILE $tmphash{$tmpkey4},"\n";
		    }
		}
	    }else{
	    }
	}
	close FILE;
	open (FILE,"$file6") or die "cannot open $file6\n";
	undef %Bfinal_chimeric_hash;
	while (<FILE>){
	    my @final_array=split/\s+/,$_;
	    $Bfinal_chimeric_hash{$final_array[0]}+=$final_array[4];
	}
	close FILE;
	my $final_key;
	open (HFILE,">$key/$tmp_name.tmpids") or die "cannot open $key/$tmp_name.tmpids";
	foreach $final_key (sort keys %Bfinal_chimeric_hash){
	    print HFILE $final_key,"\n" if $Bfinal_chimeric_hash{$final_key}>=0.8;
	}
	close HFILE;
#	$cmd="awk '{print \$1}' $file6 |awk '!a[\$1]++' >$key/$tmp_name.ids";
	$cmd="awk '!a[\$1]++' $key/$tmp_name.tmpids  >$key/$tmp_name.ids";
	system $cmd;
	push (@group1,"$tmp_name:$key/$tmp_name.ids");
    }
}
#Lable:

    if ($skip_step2 eq "Y"){
#	chdir $indenfity_doc;
	foreach $key (sort keys %focus_species_name){
	    foreach $tmp_name (sort keys %protein_path){
		push (@group1,"$tmp_name:$key/$tmp_name.ids");
		my $file4=$key."2".$tmp_name."candidate_chimerics_long_copy.ids";
		push (@group2,"$tmp_name:$key/$file4");

	    }
	}
}
#exit;
###################################################################################################################
#==# treat tree, to get the ingroup, midgroup and outgroup
# the replaced variable like "_a1, _b1" was storged in the %expand, using the $expand{"_a1"} to get it
my @iteration_level=("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z");
my %group_hash;
my %group_detail;
my @species_in_edit_tree;
&iteration_tree($phylogeny,$phylogeny,0);
unshift (@species_in_edit_tree,$phylogeny);
my %expand;

foreach my $key (sort keys %group_detail){
    &print_group_detail($group_detail{$key},$key,$key);
}

my $species_key;
my @ingroup_names;
my @midgroup_names;
my @outgroup_names;
open (ALLDETAILS,">all_details.txt") or die "wrong, cannot open all_details.txt\n";
foreach $species_key (sort keys %focus_species_name){#处理好chimeric gene在每一个物种对应的long copy和short copy之后，根据动态产生的ingroup和midgroup，对结果进行分类
    my $file7=$species_key."final_chimeric_gene_ids";
    open (FINALC,">$file7") or die "cannot open file $species_key"."final_chimeric_gene_ids";
    my $group2_key;
    my $group1_key;
    my %all_chimeric_hash;
    my %species_chimeric_hash;
    my %candidate_gene_orthologous;
    foreach $group1_key (@group1){
	my ($species,$file_name)=split/\:/,$group1_key;
	my @ttmp=split/\//,$file_name;
	next if $species_key ne $ttmp[0];#因为group2包含所有文件，所以，把不属于目前focus species的文件排除
	open (FILE,$file_name) or die "cannot open file $file_name\n\n\n";
	while ($info=<FILE>){
	    chomp $info;
	    $all_chimeric_hash{$info}++;
	    $species_chimeric_hash{$info}{$species}=1;
	}
	close FILE;
    }
    my $tmp_chimeric;
    foreach $tmp_chimeric (sort keys %all_chimeric_hash){
	foreach (@intree_check){
	    $species_chimeric_hash{$tmp_chimeric}{$_}=0 if not exists $species_chimeric_hash{$tmp_chimeric}{$_};
	    $candidate_gene_orthologous{$tmp_chimeric}{$_}=0;
	}
    }
    foreach $group2_key (@group2){
	my ($species,$file_name)=split/\:/,$group2_key;
	my @ttmp=split/\//,$file_name;
	next if $species_key ne $ttmp[0];#因为group2包含所有文件，所以，把不属于目前focus species的文件排除
	open (FILE,$file_name) or die "cannot open file $file_name\n\n\n";
	while ($info=<FILE>){
	    chomp $info;
	    @lines=split/\s+/,$info;
	    if (exists $all_chimeric_hash{$lines[0]}){
		$candidate_gene_orthologous{$lines[0]}{$species}=1;
	    }
	}
	close FILE;
    }

    my $final_rel=$species_key."_final_chimeric_genes.html";
    open (FINAL,">$final_rel") or die "wrong, cannot open $final_rel\n";
    print FINAL <<END;
    <head>
	<title>Chimeric genes at $species_key</title>
	<meta http-equiv="content-type" content="text/html; charset=UTF-8" />
	<meta name="keywords" content="chimeric gene, whole genome level scan" />
	<meta name="description" content="chimeric gene scan results" />
	<meta name="author" content="Chengjun Zhang" />
	<meta name="copyright" content="Copyright @ Long lab" />
	</head><body>
	<STYLE type="text/css">
	    td {white-space:nowrap}
        </STYLE>
END


foreach (@species_in_edit_tree){
    undef @ingroup_names;
    undef @midgroup_names;
    undef @outgroup_names;

	my @edit_tree_detail=split/\(|\,|\)|\;/,$_;
	my $tmp_name;
	foreach $tmp_name (@edit_tree_detail){
	    next if $tmp_name eq "";
	    next if exists $outgroup_species_name{$tmp_name};
	    if ($tmp_name eq $species_key){
		push (@ingroup_names,$tmp_name);
	    }elsif(substr($tmp_name,0,1) eq "_"){
		if ($expand{$tmp_name} =~/$species_key/){
		    @ingroup_names=split/\=/,$expand{$tmp_name};
		}else{
		    my @tmp_midgroup_names=split/\=/,$expand{$tmp_name};
		    @midgroup_names=(@midgroup_names,@tmp_midgroup_names);
		}
	    }else{
		@midgroup_names=(@midgroup_names,$tmp_name);
	    }
	}
    my $check_mid_is_out;
    my @new_midgroup_names;
    foreach $check_mid_is_out (@midgroup_names){
	next if exists $outgroup_species_name{$check_mid_is_out};
	push (@new_midgroup_names,$check_mid_is_out);
    }
    undef @midgroup_names;
    @midgroup_names=@new_midgroup_names;
	my $check_midgroup=join("",@midgroup_names);
	while($check_midgroup=~s/\s+//g){}
	if ($check_midgroup eq ""){
	    undef @ingroup_names;
	    undef @midgroup_names;
	    undef @outgroup_names;
	last;
	}
###################################################################################################################
#==# to findout the chimeric gene based on the phylogeny
	
	my @ingroup;
	foreach (@ingroup_names){
	    next if $_ eq "";
	    push (@ingroup,$_);
	}
	my @midgroup;
	foreach (@midgroup_names){
	    next if $_ eq "";
	    push (@midgroup,$_);
	}
	my @outgroup;
	my $tmp_outgroup;
	foreach $tmp_outgroup (sort keys %outgroup_species_name){
	    push (@outgroup,$tmp_outgroup);
	}
	my $ingroup_species=@ingroup;
	my $outgroup_species=@outgroup;
	my $midgroup_species=@midgroup;
	my $chimeric_ids;
	my $info;

	my %ingroup;
	my %ingroup_long;
	my %midgroup;
	my %midgroup_long;
	my %outgroup;
	my %outgroup_long;
	#print FINAL "INGROUP @ingroup\n";
	#print FINAL "MIDGROUP @midgroup\n";
	#print FINAL "OUTGROUP @outgroup\n";
	print ALLDETAILS "INGROUP @ingroup\n";
	print ALLDETAILS "MIDGROUP @midgroup\n";
	print ALLDETAILS "OUTGROUP @outgroup\n";
        my %short_copy;
        my %long_copy;
	print FINAL "<div><h3>$_</h3></div>\n";
	print FINAL qq!<div><table border="1" cellspacing="1" cellpadding="1" width="100%">\n!;
	print FINAL <<END;
	<TR><td>chimeric gene id</td><td colspan=$ingroup_species>INGROUP</td><td colspan=$midgroup_species>MIDGROUP</td><td colspan=$outgroup_species>OUTGROUP</td></tr>\n<TR bgcolor="#F0FFFF"><td>&nbsp;</td>
END
foreach (@ingroup){
    next if $_ eq "";
    print FINAL <<END;
    <td >$_</td>
END
}
foreach (@midgroup){
    next if $_ eq "";
    print FINAL <<END;
    <td >$_</td>
END
}
foreach (@outgroup){
    next if $_ eq "";
    print FINAL <<END;
    <td >$_</td>
END
}
	print FINAL "</tr>\n";
    foreach $key (sort keys %all_chimeric_hash){
	    $ingroup{$key}=0;
            $midgroup{$key}=0;
            $outgroup{$key}=0;
            $ingroup_long{$key}=0;
            $midgroup_long{$key}=0;
            $outgroup_long{$key}=0;
	    foreach (@ingroup){
		next if $_ eq "";
		$ingroup{$key}+=$species_chimeric_hash{$key}{$_};
		$ingroup_long{$key}+=$candidate_gene_orthologous{$key}{$_};
	    }
	    foreach (@midgroup){
		next if $_ eq "";
		$midgroup{$key}+=$species_chimeric_hash{$key}{$_};
		$midgroup_long{$key}+=$candidate_gene_orthologous{$key}{$_};
	    }
	    foreach (@outgroup){
		next if $_ eq "";
		$outgroup{$key}+=$species_chimeric_hash{$key}{$_};
		$outgroup_long{$key}+=$candidate_gene_orthologous{$key}{$_};
	    }
#	    print  "$key \t$ingroup{$key}\t $ingroup_long{$key}\t$midgroup{$key}\t$midgroup_long{$key}\t$outgroup{$key}\t$outgroup_long{$key}\n";
	    if ($ingroup{$key}>=($ingroup_species-1)){# or ($ingroup{$key} eq 1 and $ingroup_species eq 1)){
                $short_copy{$key}=1;
            }else{
                $short_copy{$key}=0;
            }
            if ($ingroup_long{$key}>=($ingroup_species-1)){# or ($ingroup_long{$key} eq 1 and $ingroup_species eq 1)){
                $long_copy{$key}=1;
            }else{
                $long_copy{$key}=0;
            }
	    if ($ingroup_species eq 1){
		$short_copy{$key}=$ingroup{$key};
		$long_copy{$key}=$ingroup_long{$key};
	    }
	    if ($midgroup_species eq 1){
		$short_copy{$key}.=$midgroup{$key};
                $long_copy{$key}.=$midgroup_long{$key};

	    }else{

		if ($midgroup{$key}>=($midgroup_species-1)){# or ($midgroup{$key} eq 1 and $midgroup_species eq 1)){
		    $short_copy{$key}.=1;
		}else{
		    $short_copy{$key}.=0;
		}
		if ($midgroup_long{$key}>=($midgroup_species-1)){# or ($midgroup_long{$key} eq 1 and $midgroup_species eq 1)){
		    $long_copy{$key}.=1;
		}else{
		    $long_copy{$key}.=0;
		}
	    }
	    if ($outgroup_species eq 1){
		$short_copy{$key}.=$outgroup{$key};
                $long_copy{$key}.=$outgroup_long{$key};
	    }else{
		if ($outgroup{$key}>=($outgroup_species-1)){# or ($outgroup{$key} eq 1 and $outgroup_species eq 1)){
		    $short_copy{$key}.=1;
		}else{
		    $short_copy{$key}.=0;
		}
		if ($outgroup_long{$key}>=($outgroup_species-1)){# or ($outgroup_long{$key} eq 1 and $outgroup_species eq 1)){
		    $long_copy{$key}.=1;
		}else{
		    $long_copy{$key}.=0;
		}
	    }
#	    print "$key\t($long_copy{$key}\t$short_copy{$key}\n";
	    if (($long_copy{$key} eq "110" and ($short_copy{$key} eq "111")) or ($long_copy{$key} eq "100" and ($short_copy{$key} eq "111" or $short_copy{$key} eq "110"))){
		print FINALC $key,"\n";
		print FINAL "<tr><td rowspan=2>$key</td>";
		foreach (@ingroup){
		    next if $_ eq "";
#		    $cmd="grep $key $species_key/$species_key"."2$_"."candidate_chimerics_long_copy.ids";
		    my $tmpfile="$species_key/$species_key"."2$_"."candidate_chimerics_long_copy.ids";
#		    print $tmpfile1,"$species_key/$species_key"."2$_"."\tWWWWWW";
#		    exit;
		    $cmd=qq!grep $key $tmpfile|awk '{if (\$1 =="$key"){print \$0}}'!;
                    my @info=`$cmd`;

		    for ($i=0;$i<=$#info;$i++){
			$info[$i]=~s/$key//;
			$info[$i].="<br>";
		    }
		    print FINAL <<END;
		    <td><a href="#">$candidate_gene_orthologous{$key}{$_}</a><br>@info</td>
END
		}
		foreach (@midgroup){
		    next if $_ eq "";
#		    $cmd="grep $key $species_key/$species_key"."2$_"."candidate_chimerics_long_copy.ids|awk '{if ($1 =="$key")}'";
		    my $tmpfile="$species_key/$species_key"."2$_"."candidate_chimerics_long_copy.ids";
		    $cmd=qq!grep $key $tmpfile|awk '{if (\$1 =="$key"){print \$0}}'!;
                    my @info=`$cmd`;
		    for ($i=0;$i<=$#info;$i++){
			$info[$i]=~s/$key//;
			$info[$i].="<br>";
		    }
		    print FINAL <<END;
		    <td><a href="#">$candidate_gene_orthologous{$key}{$_}</a><br>@info</td>
END
		}
		foreach (@outgroup){
		    next if $_ eq "";
#		    $cmd="grep $key $species_key/$species_key"."2$_"."candidate_chimerics_long_copy.ids";
		    my $tmpfile="$species_key/$species_key"."2$_"."candidate_chimerics_long_copy.ids";
		    $cmd=qq!grep $key $tmpfile|awk '{if (\$1 =="$key"){print \$0}}'!;

                    my @info=`$cmd`;
		    for ($i=0;$i<=$#info;$i++){
			$info[$i]=~s/$key//;
			$info[$i].="<br>";
		    }
		    print FINAL <<END;
		    <td><a href="#">$candidate_gene_orthologous{$key}{$_}</a><br>@info</td>
END
		}
		print FINAL "</tr><tr>";
		foreach (@ingroup){
		    next if $_ eq "";
#		    $cmd="grep $key $species_key/$species_key"."2$_"."more_than_2_short_copy.ids";
		    my $tmpfile="$species_key/$species_key"."2$_"."more_than_2_short_copy.ids";
		    $cmd=qq!grep $key $tmpfile|awk '{if (\$1 =="$key"){print \$0}}'!;

                    my @info1=`$cmd`;
		    my @info2=@info1;
		    for ($i=0;$i<=$#info1;$i++){
			$info1[$i]=~s/$key//;
			$info1[$i].="<br>";
		    }
		    print FINAL <<END;
		    <td><a href="#">$species_chimeric_hash{$key}{$_}</a><br>@info1
END
if ($species_chimeric_hash{$key}{$_} eq 1){
&draw_chimeric_region(@info1);
}else{
    my $tmpfile="$species_key/$species_key"."2$_"."candidate_chimeric_short_copy.ids";
    $cmd=qq!grep $key $tmpfile|awk '{if (\$1 =="$key"){print \$0}}'!;
    my @info1=`$cmd`;
    next unless (@info1);
    $info1[0]=~s/$key//;
    $info1[0].="<br>";
    if (@info2){
	for ($i=0;$i<=$#info2;$i++){
	    $info2[$i]=~s/$key//;
	    $info2[$i].="<br>";
	}
	&draw_chimeric_region(@info2);
    }else{
    print FINAL <<END;
    $info1[0]
END
&draw_chimeric_region($info1[0]);
    }
}
		}
		foreach (@midgroup){
		    next if $_ eq "";
		    my $tmpfile="$species_key/$species_key"."2$_"."more_than_2_short_copy.ids";
		    $cmd=qq!grep $key $tmpfile|awk '{if (\$1 =="$key"){print \$0}}'!;
                    my @info1=`$cmd`;
		    my @info2=@info1;
		    for ($i=0;$i<=$#info1;$i++){
			$info1[$i]=~s/$key//;
			$info1[$i].="<br>";
		    }
		    print FINAL <<END;
		    <td><a href="#">$species_chimeric_hash{$key}{$_}</a><br>@info1
END
if ($species_chimeric_hash{$key}{$_} eq 1){
&draw_chimeric_region(@info1);
}else{
    my $tmpfile="$species_key/$species_key"."2$_"."candidate_chimeric_short_copy.ids";
    $cmd=qq!grep $key $tmpfile|awk '{if (\$1 =="$key"){print \$0}}'!;
    my @info1=`$cmd`;
    next unless (@info1);
    $info1[0]=~s/$key//;
    $info1[0].="<br>";
    if (@info2){
	for ($i=0;$i<=$#info2;$i++){
	    $info2[$i]=~s/$key//;
	    $info2[$i].="<br>";
	}
	&draw_chimeric_region(@info2);
    }else{
    print FINAL <<END;
    $info1[0]
END
&draw_chimeric_region($info1[0]);
    }
}
		}
		foreach (@outgroup){
		    next if $_ eq "";
		    my $tmpfile="$species_key/$species_key"."2$_"."more_than_2_short_copy.ids";
		    $cmd=qq!grep $key $tmpfile|awk '{if (\$1 =="$key"){print \$0}}'!;
                    my @info1=`$cmd`;
		    my @info2=@info1;
		    for ($i=0;$i<=$#info1;$i++){
			$info1[$i]=~s/$key//;
			$info1[$i].="<br>";
		    }
		    print FINAL <<END;
		    <td><a href="#">$species_chimeric_hash{$key}{$_}</a><br>@info1
END
if ($species_chimeric_hash{$key}{$_} eq 1){
&draw_chimeric_region(@info1);
}else{
    my $tmpfile="$species_key/$species_key"."2$_"."candidate_chimeric_short_copy.ids";
    $cmd=qq!grep $key $tmpfile|awk '{if (\$1 =="$key"){print \$0}}'!;
    my @info1=`$cmd`;
    next unless (@info1);
    $info1[0]=~s/$key//;
    $info1[0].="<br>";
    if (@info2){
	for ($i=0;$i<=$#info2;$i++){
	    $info2[$i]=~s/$key//;
	    $info2[$i].="<br>";
	}
	&draw_chimeric_region(@info2);
    }else{
    print FINAL <<END;
    $info1[0]
END
&draw_chimeric_region($info1[0]);
    }
}
		}
		print FINAL "</tr>";
	    }

	    
	}
    print FINAL "</table></div>";
    
}
    close FINAL;
    close FINALC;
}

sub draw_chimeric_region(){
    my @draw=@_;
    my @tmpdraw;
    my @finaldraw;
    my $percent;
#    print "WW @draw\n";
    for ($i=0;$i<=$#draw;$i++){
	@tmpdraw=split/\s+/,$draw[$i];
	push (@finaldraw,$tmpdraw[2]);
	push (@finaldraw,$tmpdraw[3]);
#	print "DD $tmpdraw[3] $tmpdraw[2]\n";
    }
    undef @tmpdraw;
    @tmpdraw=sort @finaldraw;
    print FINAL qq!<table border="0" cellspacing="0" cellpadding="0" width="100%"><tr bgcolor="red">!;
    for ($i=0;$i<=$#tmpdraw;$i++){
	if ($i eq "0"){
	    $percent=int($tmpdraw[$i]*100-1);
	}else{
	    $percent=int(($tmpdraw[$i]-$tmpdraw[$i-1])*100);
	}
	if ($i/2 eq int ($i/2)){
	    print FINAL qq!<td width="$percent%">&nbsp;</td>!;
	}else{
	    print FINAL qq!<td width="$percent%" bgcolor="green">&nbsp;</td>!;
	}
    }
#   print "EE @tmpdraw $tmpdraw[$#tmpdraw]\n";
    
    $percent=100-int($tmpdraw[$#tmpdraw]*100);
    print FINAL qq!<td width="$percent%">&nbsp;</td>!;
    print FINAL "</tr></table></td>";
#    exit;
}

sub iteration_tree(){
    my ($tree,$original_tree,$iteration_times)=@_;
    while ($tree=~s/\(\(/\(/){}
    while ($tree=~s/\)\)/\)/){}
    my $i=0;
    while ($tree=~s/\(\w+,\w+\)//){
        $i++;
        $group_hash{$&}="_".$iteration_level[$iteration_times].$i;
        $group_detail{"_".$iteration_level[$iteration_times].$i}=$&;
    }
    $iteration_times++;
    my $key2;
    foreach $key2 (sort keys %group_hash){
        $original_tree=~s/\($key2\)/$group_hash{$key2}/;
    }
    push (@species_in_edit_tree,$original_tree);
    &iteration_tree($original_tree,$original_tree,$iteration_times) if $original_tree=~/\(|\)/;
}


sub print_group_detail(){
    my ($detail,$key,$key1)=@_;
    my @info=split/\(|\,|\)/,$detail;
    foreach (@info){
        while ($_=~s/\s+//g){}
        next if $_ eq "";
        if (substr($_,0,1) eq "_"){
            &print_group_detail($group_detail{$_},$_,$key1);
        }else{
            $expand{$key1}.="=".$_;
        }
    }
}


sub tree(){
    my ($tree11)=@_;
    my $tmptree=$tree11;
    $tmptree=~s/\(|\)|\;//gi;
    my @intrees=split/\,/,$tmptree;
    return @intrees;
}

sub read_fasta_file_new()
{
    my %hash;
    my ($file_name,$sep)=@_;
    my $aa;
    my $id;
    my @lines;
    open (PEP,"$file_name") or die "wrong";
    my $first_id=<PEP>;
    chomp $first_id;
    $first_id=~s/>//;
    if ($sep eq "dot"){
        @lines=split/\./,$first_id;
    }elsif($sep eq "shu"){
        @lines=split/\|/,$first_id;
    }elsif($sep eq "space"){
        @lines=split/\s/,$first_id;
    }elsif($sep eq "join"){
        @lines=split/\s+|\t/,$first_id;
        my $tmp=join("-",@lines);
        $lines[0]=$tmp;
    }else{
        @lines=($first_id);
    }

    $id=$lines[0];
    while (<PEP>){
        chomp;
        if ($_=~s/>//){
            $aa=~s/\*//;
            $hash{$id}=$aa;
            $aa="";
            if ($sep eq "dot"){
                @lines=split/\./,$_;
            }elsif($sep eq "shu"){
                @lines=split/\|/,$_;
            }elsif($sep eq "space"){
                @lines=split/\s/,$_;
            }elsif($sep eq "join"){
                @lines=split/\s+|\t/,$_;
                my $tmp=join("-",@lines);
                $lines[0]=$tmp;
            }else{
                @lines=($_);
            }
            $id=$lines[0];
        }else{
            $aa.=$_;
        }
    }

    $aa=~s/\*//;
    $hash{$id}=$aa;
    return %hash;
}

sub scan_chimeric()
{
#是file3的第一列,file3是存储着的blastp结果======似乎这里会混淆？？
#这个子程序，主要是生成了一个long-copy id文件，以及一个short-copy id文件，并没有保存任何哈希或者数组继续用于程序
#这个子程序的%hash_scan数据继承自主程序，其包含三个key，第一个为关心的chimeric gene id，第二个为对应的short copy id,第三个为hit对应于chimeric gene id的起止和对应于short copy id的起止

    my ($tmpid)=@_;
    my %iteration;
    my %size;
    undef %size;
    undef %iteration;
    my $tmpkey;
    
    foreach $tmpkey (sort keys %{$hash_scan{$tmpid}}){#tmpid是chimeric gene id,$tmpkey是对应的short copy的id
	print $tmpkey if not exists  $query_aa{$tmpid};
	if (length($query_aa{$tmpid}) ne 0 and length($target_aa{$tmpkey}) ne 0){
	    $qstart_locat=int($hash_scan{$tmpid}{$tmpkey}{q_start}/length($query_aa{$tmpid})*1000+0.5)/1000;
	    $qstop_locat=int($hash_scan{$tmpid}{$tmpkey}{q_stop}/length($query_aa{$tmpid})*1000+0.5)/1000;
	    $tstart_locat=int($hash_scan{$tmpid}{$tmpkey}{t_start}/length($target_aa{$tmpkey})*1000+0.5)/1000;
	    $tstop_locat=int($hash_scan{$tmpid}{$tmpkey}{t_stop}/length($target_aa{$tmpkey})*1000+0.5)/1000;
	    $final_score=$qstop_locat- $qstart_locat+0.001;#+$tstop_locat-$tstart_locat+0.001;
	}
	if ($final_score<0.8){#原先为相加小于1.6，后来觉得实在不如直接小于0.8来得直接
	    $iteration{$tmpkey}{start}=$qstart_locat;
	    $iteration{$tmpkey}{stop}=$qstop_locat;
	    $iteration{$tmpkey}{tstart}=$tstart_locat;
	    $iteration{$tmpkey}{tstop}=$tstop_locat;
	    $size{$tmpkey}=$qstop_locat-$qstart_locat;
	}else{
            print COPIES "$tmpid\t$tmpkey\t$qstart_locat\t$qstop_locat\t",$qstop_locat-$qstart_locat,"\t$tstart_locat\t$tstop_locat\n";
	}
    }

    my $tmpkey2;
    foreach $tmpkey (sort {$size{$b}<=>$size{$a}} keys %size){#这个key是short copy id
	next if not exists $size{$tmpkey};
	foreach $tmpkey2 (sort {$size{$b}<=>$size{$a}} keys %size){
#在这里，也会涉及点的问题
#这里的点，是short copy的id，似乎没有太多影响
	    my @ttmp1=split/\./,$tmpkey;
	    my @ttmp2=split/\./,$tmpkey2;
	    next if $ttmp1[0] eq $ttmp2[0];
	    if ($iteration{$tmpkey}{start} and $iteration{$tmpkey}{stop} and $iteration{$tmpkey2}{start} and $iteration{$tmpkey2}{stop}){
                $overlap=&compare_cover($iteration{$tmpkey}{start}, $iteration{$tmpkey}{stop},$iteration{$tmpkey2}{start}, $iteration{$tmpkey2}{stop});
            }
            if( $overlap eq 0){
            }else{
                delete $size{$tmpkey2};
                delete $iteration{$tmpkey2};
            }
        }
    }
    foreach $tmpkey (sort {$size{$b}<=>$size{$a}} keys %size){
        print SCAN "$tmpid\t$tmpkey\t$iteration{$tmpkey}{start}\t$iteration{$tmpkey}{stop}\t$size{$tmpkey}\t$iteration{$tmpkey}{tstart}\t$iteration{$tmpkey}{tstop}\n" if $size{$tmpkey}>=0.2;#仅留下长度大约chimeric gene 0.2以上的copy，用于后续分析
    }
}
sub compare_cover(){#用于比较4个参数，a,b;c,d是否存在overlap，当overlap大于等于20bp的时候，就可以认为该hit包含了目标序列——intron丢失
    my ($a,$b,$c,$d)=@_;
    ($a,$b)=&switch($a,$b);
    ($c,$d)=&switch($c,$d);
    $overlap=0;
    if (($c>=$a and $c<=$b) and $d<=$b){
        $overlap=$d-$c;
    }elsif (($c>=$a and $c<=$b) and $d>$b){
        $overlap=$b-$c;
    }elsif($d>=$a and $d<=$b){
        $overlap=$d-$a;
    }elsif ($a>$c and $b<$d){
        $overlap=$b-$a;
    }else{
        $overlap=0;
    }
    return $overlap;
}
sub switch(){
    my ($a,$b)=@_;
    if ($a>$b){
        $t=$a;
        $a=$b;
        $b=$t;
    }
    @new=($a,$b);
    return (@new);
}


