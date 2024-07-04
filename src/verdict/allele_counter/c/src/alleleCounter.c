/**   LICENSE
* Copyright (c) 2014-2020 Genome Research Ltd.
*
* Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
*
* This file is part of alleleCount.
*
* alleleCount is free software: you can redistribute it and/or modify it under
* the terms of the GNU Affero General Public License as published by the Free
* Software Foundation; either version 3 of the License, or (at your option) any
* later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
* details.
*
* You should have received a copy of the GNU Affero General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <bam_access.h>
#include "dbg.h"

static int min_base_q = 20;
static int min_map_q = 35;
static char *hts_file;
static char *loci_file;
static char *out_file;
static char *ref_file;
static char *contig = NULL;
static int inc_flag = 3; //Paired, proper pair
static int exc_flag = 3852; // Read unmapped, Mate unmapped, Secondary alignment, Fails QC, Duplicate, Supplementary alignment
static int snp6 = 0;
static int is_10x = 0;
static int is_dense = 0;

int check_exist(char *fname){
	FILE *fp;
	if((fp = fopen(fname,"r"))){
		fclose(fp);
		return 1;
	}
	return 0;
}

void alleleCounter_print_usage (int exit_code){

	printf ("Usage: alleleCounter -l loci_file.txt -b sample.bam -o output.txt [-m int] [-r ref.fa.fai]\n\n");
  printf (" -l  --loci-file [file]           Path to loci file.\n");
  printf (" -b  --hts-file [file]            Path to sample HTS file.\n");
  printf (" -o  --output-file [file]         Path write output file.\n\n");

	printf ("Optional\n");
	printf (" -r  --ref-file [file]           Path to reference fasta index file.\n");
	printf ("                                 NB. If cram format is supplied via -b and the reference listed in the cram header\n");
	printf ("                                     can't be found alleleCounter may fail to work correctly.\n");
	printf (" -m  --min-base-qual [int]       Minimum base quality [Default: %d].\n",min_base_q);
	printf (" -q  --min-map-qual [int]        Minimum mapping quality [Default: %d].\n",min_map_q);
	printf (" -c  --contig [string]           Limit calling to named contig.\n");
	printf (" -d  --dense-snps                Improves performance where many positions are close together \n");
	printf (" -x  --is-10x                    Enables 10X processing mode.\n");
	printf ("                                   In this mode the HTS input file must be a cellranger produced BAM file.  Allele\n");
	printf ("                                   counts are then given on a per-cellular barcode basis, with each count representing\n");
	printf ("                                   the consensus base for that UMI. \n");
	printf ("                                 by iterating through bam file rather than using a 'fetch' approach.\n");
	printf (" -f  --required-flag [int]       Flag value of reads to retain in allele counting default: [%i].\n",inc_flag);
  printf ("                                 N.B. if the proper-pair flag is are selected, alleleCounter will assume paired-end\n");
  printf ("                                 and filter out any proper-pair flagged reads not in F/R orientation.");
	printf (" -F  --filtered-flag [int]       Flag value of reads to exclude in allele counting default: [%i].\n",exc_flag);
	printf (" -v  --version                   Display version number.\n");
	printf (" -h  --help                      Display this usage information.\n\n");
  exit(exit_code);
}

void alleleCounter_print_version (int exit_code){
  printf ("%s\n",ALLELECOUNTER_VERSION);
  exit(exit_code);
}

void alleleCounter_setup_options(int argc, char *argv[]){
  ref_file = NULL;
	const struct option long_opts[] =
	{
             	{"loci-file", required_argument, 0, 'l'},
             	{"hts-file", required_argument, 0, 'b'},
             	{"ref-file", required_argument, 0, 'r'},
             	{"output-file",required_argument , 0, 'o'},
             	{"min-base-qual", required_argument, 0, 'm'},
							{"min-map-qual", required_argument, 0, 'q'},
							{"is-snp6", required_argument, 0, 's'},
							{"is-10x", required_argument, 0, 'x'},
							{"contig", required_argument, 0, 'c'},
							{"dense-snps", no_argument, 0, 'd'},
							{"required-flag", required_argument, 0, 'f'},
							{"filtered-flag", required_argument, 0, 'F'},
							{"version", no_argument, 0, 'v'},
             	{"help", no_argument, 0, 'h'},
             	{ NULL, 0, NULL, 0}
   }; //End of declaring opts

   int index = 0;
   int iarg = 0;

   //Iterate through options
   while((iarg = getopt_long(argc, argv, "f:F:l:b:m:o:q:r:c:hdsvx", long_opts, &index)) != -1){
   	switch(iarg){
   		  case 'h':
         	alleleCounter_print_usage(0);
         	break;

        case 'v':
          alleleCounter_print_version(0);
          break;

      	case 'l':
      		loci_file = optarg;
      		break;

      	case 'm':
      		min_base_q = atoi(optarg);
      		break;

        case 'r':
          ref_file = optarg;
          break;

      	case 'q':
      		min_map_q = atoi(optarg);
      		break;

      	case 'b':
      		hts_file = optarg;
      		break;

      	case 'o':
      		out_file = optarg;
      		break;

      	case 's':
      		snp6 = 1;
      		break;

         case 'x':
            is_10x = 1;
            break;

      	case 'c':
      	  contig = optarg;
          break;

        case 'd':
          is_dense = 1;
          break;

        case 'f':
          inc_flag = atoi(optarg);
          break;

        case 'F':
          exc_flag = atoi(optarg);
          break;

				case '?':
          alleleCounter_print_usage (1);
          break;

      	default:
      		alleleCounter_print_usage (1);

   	}; // End of args switch statement

   }//End of iteration through options

   if(check_exist(loci_file) != 1){
   	printf("Loci file %s does not appear to exist.\n",loci_file);
   	alleleCounter_print_usage(1);
   }
   if(check_exist(hts_file) != 1){
   	printf("HTS file %s does not appear to exist.\n",hts_file);
   	alleleCounter_print_usage(1);
   }
   if(ref_file){
    if( check_exist(ref_file) != 1){
   	  printf("Reference file provided %s does not appear to exist.\n",ref_file);
   	  alleleCounter_print_usage(1);
   	 }
   }
   return;
}

int print_loci_head(FILE *output){
	int chk = fprintf(output,"#CHR\tPOS\tCount_A\tCount_C\tCount_G\tCount_T\tGood_depth\n");
	return chk;
}

int print_snp6_header(FILE *output){
	int chk = fprintf(output,"#CHR\tPOS\tCount_Allele_A\tCount_Allele_B\tGood_depth\n");
	return chk;
}

int print_10x_header(FILE *output){
	int chk = fprintf(output,"#CHR\tPOS\tBarcode\tCount_A\tCount_C\tCount_G\tCount_T\tGood_depth\n");
   return chk;
}

int print_header(FILE *output, int snp6){
	if(snp6 == 1){
		return print_snp6_header(output);
	}else if(is_10x == 1){
     return print_10x_header(output);
   }else{
		return print_loci_head(output);
	}
}

int calculateAlleleCount(int a_cnt, int c_cnt, int g_cnt, int t_cnt, int allele){
	switch (allele){
		case 'A':
			return a_cnt;

		case 'C':
			return c_cnt;

		case 'G':
			return g_cnt;

		case 'T':
			return t_cnt;

		default:
			return -1;
	};
}

int print_loci_section(FILE *output, char *chr, int pos, int a_cnt, int c_cnt, int g_cnt, int t_cnt, int depth){
	assert(output !=NULL);
	return (fprintf(output,"%s\t%d\t%d\t%d\t%d\t%d\t%d\n",chr,pos,a_cnt,c_cnt,g_cnt,t_cnt,depth));
}

int print_snp6_section(FILE *output, char *chr, int pos, int allele_a, int allele_b, int depth){
	assert(output !=NULL);
	return (fprintf(output,"%s\t%d\t%d\t%d\t%d\n",chr,pos,allele_a,allele_b,depth));
}

int print_section(FILE *output, char *chr, int pos, int a_cnt, int c_cnt, int g_cnt,
															int t_cnt, int depth, int snp6, char allele_A, char allele_B){
	if(snp6 == 1){
		int all_a_cnt = calculateAlleleCount(a_cnt, c_cnt, g_cnt, t_cnt, allele_A);
		check(all_a_cnt>=0,"Error getting A Allele count '%c'",allele_A);
		int all_b_cnt = calculateAlleleCount(a_cnt, c_cnt, g_cnt, t_cnt, allele_B);
		check(all_b_cnt>=0,"Error getting B Allele count '%c'",allele_B);
		return print_snp6_section(output, chr, pos, all_a_cnt, all_b_cnt, depth);
	}else{
		return print_loci_section(output, chr, pos, a_cnt, c_cnt, g_cnt,t_cnt, depth);
	}
	error:
		return -1;
}

int get_position_info_from_file(char *line, loci_stats *stats, int snp6, int i){
	int chr_d = 0;

	if(snp6==1){
		int chk = sscanf(line,"%d%*[ \t]%d%*[ \t]%*s%*[ \t]%*s%*[ \t]%c%*[ \t]%c",&chr_d,&(stats->pos),&(stats->allele_A),&(stats->allele_B));
		if(chk == 2){
			int try = sprintf(stats->chr,"%d",chr_d);
			check(try >0,"Error trying to convert chromosome name '%d'to string.",chr_d);
		}else{
			//Try again but a string match
			chk = sscanf(line,"%s%*[ \t]%d%*[ \t]%*s%*[ \t]%*s%*[ \t]%c%*[ \t]%c",stats->chr,&(stats->pos),&(stats->allele_A),&(stats->allele_B));
			check(chk==4,"Error attempting string match of allele position info from SNP6 line %s.",line);
		}
		check(chk==2,"Error parsing SNP6 file line number %d: '%s'.",i,line);
	}else{
		int chk = sscanf(line,"%d%*[ \t]%d",&chr_d,&(stats->pos));
		if(chk == 2){
			int try = sprintf(stats->chr,"%d",chr_d);
			check(try >0,"Error trying to convert chromosome name '%d'to string.",chr_d);
		}else{
			//Try again but a string match
			chk = sscanf(line,"%s%*[ \t]%d",stats->chr,&(stats->pos));
			check(chk==2,"Error parsing loci file line number %d as a string match: '%s'.",i,line);
		}
		check(chk==2,"Error parsing loci file line number %d: '%s'.",i,line);
	}
	return 0;
	error:
		return -1;
}

int line_count (char *file_path){
  FILE *f = fopen(file_path,"r");
  int line_count = 0;
  check(f != NULL, "Error opening file '%s' to count lines.",file_path);
  char rd[ 5000 ];
	while(fgets(rd, sizeof(rd), f) != NULL){
    line_count++;
  }
  fclose(f);
  return line_count;
	error:
  if(f) fclose(f);
  return -1;
}

int sort_loci_stats(const void *a1, const void *b1){
	loci_stats *a = *(loci_stats * const *)a1;
	loci_stats *b = *(loci_stats * const *)b1;
	int res = strcmp(a->chr,b->chr);
	if(res==0){
		if(a->pos == b->pos){
			return 0;
		}else{
			return a->pos < b->pos ? -1 : 1;
		}
	}else{
		return res;
	}
}

int init_base_counts(loci_stats *stats){
	stats->base_counts = malloc(sizeof(int) * 4);
	check_mem(stats->base_counts);
	stats->base_counts[0] = 0;
	stats->base_counts[1] = 0;
	stats->base_counts[2] = 0;
	stats->base_counts[3] = 0;
	return 0;
	error:
	return 1;
}

loci_stats ** read_locis_from_file(char *loci_file, int *line_cnt){
	FILE *loci_in = NULL;
	loci_stats **stats= NULL;
	*line_cnt = line_count(loci_file);
	check(*line_cnt>=0,"Error counting lines in loci file: %s",loci_file);

	stats = malloc(sizeof(loci_stats*)*(*line_cnt));
	check_mem(stats);
	//Open loci file
  loci_in = fopen(loci_file,"r");
  check(loci_in != NULL, "Error opening loci file %s for reading.",loci_file);
	int i=0;
	char line[2048];
	while ( fgets(line,sizeof(line),loci_in) != NULL ){
		stats[i] = malloc(sizeof(loci_stats));
		check_mem(stats[i]);
		stats[i]->chr = malloc(sizeof(char)*2048);
		check_mem(stats[i]->chr);
		stats[i]->base_counts = NULL;
		int check = get_position_info_from_file(line,stats[i],snp6,i);
		check(check==0,"Error trying to fetch position from file at line %d.",i);
		check = init_base_counts(stats[i]);
		check(check==0,"Error initialising base counts %d.",i);
		i++;
	}
	int size = *line_cnt;
	qsort(stats,size,sizeof(loci_stats*),&sort_loci_stats);
	fclose(loci_in);
	return stats;
	error:
		if(stats) {
			int j=0;
			for(j=0;j<*line_cnt;j++){
				if(stats[j]){
					free(stats[j]->chr);
					free(stats[j]);
				}
			}
			free(stats);
		}
		if(loci_in) fclose(loci_in);
		return NULL;
}

int main(int argc, char *argv[]){
	loci_stats **locis = NULL;
	//Get the options commandline
	alleleCounter_setup_options(argc,argv);
	//Set the min base and mapping quality.
	bam_access_min_base_qual(min_base_q);

	bam_access_min_map_qual(min_map_q);

	bam_access_inc_flag(inc_flag);

	bam_access_exc_flag(exc_flag);

	//Open output file for writing
	FILE *output = fopen(out_file,"w");
  check(output != NULL, "Error opening file %s for write.",out_file);
	int chk = print_header(output,snp6);
	check(chk >= 0,"Error trying to write header '%s'.",out_file);
	//Open bam file and iterate through chunks until we reach the cutoff.
	chk = -1;

	chk = bam_access_openhts(hts_file,ref_file);
	check(chk == 0,"Error trying to open sequence/index files '%s'.",hts_file);
	int loci_count=0;
	fprintf(stderr,"Reading locis\n");
	locis = read_locis_from_file(loci_file,&loci_count);
	fprintf(stderr,"Done reading locis\n");

	check(locis!=NULL,"Error reading loci_stats from file.");
   if(is_10x){
     fprintf(stderr,"Using 10X processing mode.\n");
     if(!is_dense){
		  int j=0;
		  for(j=0;j<loci_count;j++){
		  	int ret = bam_access_get_position_base_counts(locis[j]->chr,locis[j]->pos,locis[j],is_10x,output);
		  	check(ret==0,"Error retrieving stats from bam file for position %s:%d",locis[j]->chr,locis[j]->pos);
		  	free(locis[j]->chr);
		  	if(locis[j]->base_counts) free(locis[j]->base_counts);
		  	free(locis[j]);
        }
     }else{
  		 fprintf(stderr,"Multi pos start:\n");
  		 int ret = bam_access_get_multi_position_base_counts(locis, loci_count,is_10x,output);
  		 check(ret==0,"Error scanning through bam file for loci list with dense snps.");
     }
   }else{
    if(is_dense){
  		fprintf(stderr,"Multi pos start:\n");
  		int ret = bam_access_get_multi_position_base_counts(locis, loci_count,is_10x,output);
  		check(ret==0,"Error scanning through bam file for loci list with dense snps.");
  		int j=0;
  		for(j=0;j<loci_count;j++){
  			int depth = locis[j]->base_counts[0]+locis[j]->base_counts[1]+locis[j]->base_counts[2]+locis[j]->base_counts[3];
        int check_print = print_section(output,locis[j]->chr,locis[j]->pos,locis[j]->base_counts[0],
                    locis[j]->base_counts[1],locis[j]->base_counts[2],locis[j]->base_counts[3],depth,
                    snp6,locis[j]->allele_A,locis[j]->allele_B);
        check(check_print>0,"Error printing line to output file: %s: %d.",locis[j]->chr,locis[j]->pos);
  			free(locis[j]->chr);
  			if(locis[j]->base_counts) free(locis[j]->base_counts);
  			free(locis[j]);
  		}
    }else{
  		int j=0;
  		for(j=0;j<loci_count;j++){
  			int ret = bam_access_get_position_base_counts(locis[j]->chr,locis[j]->pos,locis[j],is_10x,output);
  			check(ret==0,"Error retrieving stats from bam file for position %s:%d",locis[j]->chr,locis[j]->pos);
        int depth = locis[j]->base_counts[0]+locis[j]->base_counts[1]+locis[j]->base_counts[2]+locis[j]->base_counts[3];
        int check_print = print_section(output,locis[j]->chr,locis[j]->pos,locis[j]->base_counts[0],
                    locis[j]->base_counts[1],locis[j]->base_counts[2],locis[j]->base_counts[3],depth,
                    snp6,locis[j]->allele_A,locis[j]->allele_B);
        check(check_print>0,"Error printing line to output file: %s: %d.",locis[j]->chr,locis[j]->pos);

  			free(locis[j]->chr);
  			if(locis[j]->base_counts) free(locis[j]->base_counts);
  			free(locis[j]);
  		}
  		free(locis);
  	 }
   }

	//Close files.
	//fclose(loci_in);
	bam_access_closehts();
	fclose(output);
	return 0;

error:
	bam_access_closehts();
	if(locis){
		int j=0;
		for(j=0;j<loci_count;j++){
				if(locis[j]){
					free(locis[j]->chr);
					if(locis[j]->base_counts) free(locis[j]->base_counts);
					free(locis[j]);
				}
		}
		free(locis);
	}
	if(output) fclose(output);
	if(hts_file) free(hts_file);
	if(out_file) free(out_file);
	if(loci_file) free(loci_file);
	return 1;
}
