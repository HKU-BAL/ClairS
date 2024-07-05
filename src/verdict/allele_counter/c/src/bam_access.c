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

#include <bam_access.h>
#include <stdlib.h>
#include <stdio.h>
#include <dbg.h>
#include <assert.h>
#include <limits.h>
#include <htslib/cram.h>
#include "khash.h"

#define PO10_LIMIT (INT_MAX/10)
KHASH_MAP_INIT_STR(strh,uint8_t)

file_holder *fholder;
int counter = -1;
int include_sw = 0;
int include_dup = 0;
int include_se = 0;
int min_base_qual = 20;
int min_map_qual = 35;
int inc_flag = 3;
int exc_flag = 3852;
int maxitercnt = 1000000000; //Overrride internal maxcnt for iterator!
//Make sure this isn't too close to the integer overflow boundary
//int maxitercnt = 100000;

typedef struct {

} plp_aux_t;

int print_10x_section(FILE *output, char *chr, int pos, int a_cnt, int c_cnt, int g_cnt, int t_cnt, int depth,char *barcode){
	assert(output !=NULL);
	return (fprintf(output,"%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n",chr,pos,barcode,a_cnt,c_cnt,g_cnt,t_cnt,depth));
}

int bam_access_openhts(char *hts_file, char *ref_file){
	assert(hts_file != NULL);
	//Assign memory for the file name etc holding struct
	fholder = malloc(sizeof(file_holder));
	check_mem(fholder);
	//Beginning and end of tmp struct for bam access
	fholder->beg = 0; fholder->end = 0x7fffffff;  // The max 32 bit integer.
	//Open a file for read from compressed bam.
	fholder->in = hts_open(hts_file, "r");
	check(fholder->in != 0,"HTS file %s failed to open.",hts_file);
  fholder->idx = sam_index_load(fholder->in,hts_file);
	check(fholder->idx != 0,"HTS index file for %s failed to open.",hts_file);
	if(ref_file){
	  int chk = hts_set_fai_filename(fholder->in, ref_file);
		check(chk==0,"Error setting fai filename %s.",ref_file);
	}else{
	  if(fholder->in->format.format == cram) log_warn("No reference file provided for a cram input file, if the reference described in the cram header can't be located this script may fail.");
	}
  //Check for generic header read method.
  fholder->head = sam_hdr_read(fholder->in);
	return 0;
error:
  if(fholder->idx) hts_idx_destroy(fholder->idx);
	if(fholder->in) hts_close(fholder->in);
	if(fholder->head) bam_hdr_destroy(fholder->head);
	if(fholder) free(fholder);
	return -1;
}

void bam_access_closehts(){
	if(fholder && fholder->idx) hts_idx_destroy(fholder->idx);
	if(fholder && fholder->in) hts_close(fholder->in);
	if(fholder && fholder->head) bam_hdr_destroy(fholder->head);
	if(fholder) free(fholder);
	return;
}

int no_of_digits(int i){
	int n,po10;

  if (i < 0) i = -i;
  n=1;
  po10=10;
  while(i>=po10)
  {
    n++;
    if (po10 > PO10_LIMIT) break;
    po10*=10;
  }
  return n;
}

//callback for bam_plp_init
static int pileup_func(void *data, bam1_t *b){
  return 0;
}

void pileupCounts(const bam_pileup1_t *pil, int n_plp, loci_stats *stats){
	khash_t(strh) *h;
	khiter_t k;
	h = kh_init(strh);
	int i=0;
	for(i=0;i<n_plp;i++){
		const bam_pileup1_t *p = pil + i;
		int qual = bam_get_qual(p->b)[p->qpos];
		uint8_t c = bam_seqi(bam_get_seq(p->b), p->qpos);
		int absent;
    k = kh_put(strh, h, bam_get_qname(p->b), &absent);
		uint8_t pre_b;
		if(!absent){ //Read already processed to get base processed (we only increment if base is different between overlapping read pairs)
			k = kh_get(strh, h, bam_get_qname(p->b));
			pre_b = kh_val(h,k);
		}else{
			//Add the value to the hash
			kh_value(h, k) = c;
		}
		if(!(p->is_del) &&  qual >= min_base_qual && (absent || pre_b != c)){
			//&& (c == 1 /*A*/|| c == 2 /*C*/|| c == 4 /*G*/|| c == 8 /*T*/)){
			//Now we add a new read pos struct to the list since the read is valid.
			//char cbase = toupper(bam_nt16_rev_table[c]);
			switch(c){
				case 1:
				stats->base_counts[0]++;
				break;

				case 2:
				stats->base_counts[1]++;
				break;

				case 4:
				stats->base_counts[2]++;
				break;

				case 8:
				stats->base_counts[3]++;
				break;

				default:
				break;

			}; // End of args switch statement */
		}
	}
	kh_destroy(strh, h);
	return;
}

//Sort pileup
int readCompare(const void *r1,const void *r2){
  char *barcode1,*barcode2;
  char *umi1,*umi2;
  int bcomp;
  bam_pileup1_t *e1 = (bam_pileup1_t *)r1;
  bam_pileup1_t *e2 = (bam_pileup1_t *)r2;
  barcode1 = bam_aux2Z(bam_aux_get(e1->b,"CB"));
  barcode2 = bam_aux2Z(bam_aux_get(e2->b,"CB"));
  bcomp = strcmp(barcode1,barcode2);
  //Is the barcode equal?
  if(bcomp==0){
    //If it is, sort on UMI
    umi1 = bam_aux2Z(bam_aux_get(e1->b,"UB"));
    umi2 = bam_aux2Z(bam_aux_get(e2->b,"UB"));
    return strcmp(umi1,umi2);
  }
  return bcomp;
}

void pileupCounts10x(const bam_pileup1_t *pil, int n_plp, loci_stats *stats,FILE *output){
  int i,j,k;
  char *barcode = NULL;
  char *umi = NULL;
  char *curr_barcode = NULL;
  char *curr_umi = NULL;
  int cnts[4] = {0};
  int cellCnts[4] = {0};
  bam_pileup1_t *p;
  int qual;
  int max_obs_reads;
  uint8_t c;
  //Make a non-constant pointer to reads
  p = (bam_pileup1_t *)pil;
  //Sort the input array
  qsort(p,n_plp,sizeof(bam_pileup1_t),readCompare);
  //Loop over sorted reads
  printf("Performing pileup of %d reads at %s %d\n",n_plp,stats->chr,stats->pos);
  for(i=0;i<n_plp;i++){
    qual = bam_get_qual(p->b)[p->qpos];
    c = bam_seqi(bam_get_seq(p->b), p->qpos);
    //Get the tags
    barcode = bam_aux2Z(bam_aux_get(p->b,"CB"));
    umi = bam_aux2Z(bam_aux_get(p->b,"UB"));
    //printf("CB=%s, UB=%s\n",barcode,umi);
    //Skip this read?
    if((p->is_del) || qual < min_base_qual){
      //printf("Skipping with %d and qual %d.\n",p->is_del,qual);
      p++;
      continue;
    }
    //First time, so we need to initialise current barcode/umi
    if(curr_umi==NULL && curr_barcode==NULL){
      curr_barcode = barcode;
      curr_umi = umi;
    }
    //Count them now we can assume they're sorted
    //Check if the UMI has changed?
    if(strcmp(umi,curr_umi)!=0){
      //Get the consensus read
      max_obs_reads = -1;
      j=-1;
      for(k=0;k<4;k++){
        if(cnts[k] == max_obs_reads){
          j=-1;
        }
        if(cnts[k] > max_obs_reads){
          max_obs_reads = cnts[k];
          j=k;
        }
      }
      //Add it to the cell level counter
      if(j<0){
         //printf("No consensus allele: %d,%d,%d,%d\n",cnts[0],cnts[1],cnts[2],cnts[3]);
      }else{
        cellCnts[j]++;
      }
      //Re-zero and store new current UMI
      cnts[0]=cnts[1]=cnts[2]=cnts[3]=0;
      curr_umi = umi;
    }
    //Has the barcode changed?
    if(strcmp(barcode,curr_barcode)!=0){
      //The barcode has changed, so print the old one (assuming we found something to use)
      if(cellCnts[0]+cellCnts[1]+cellCnts[2]+cellCnts[3]>0)
        print_10x_section(output,stats->chr,stats->pos,cellCnts[0],cellCnts[1],cellCnts[2],cellCnts[3],cellCnts[0]+cellCnts[1]+cellCnts[2]+cellCnts[3],curr_barcode);
      //Re-zero counters
      cellCnts[0]=cellCnts[1]=cellCnts[2]=cellCnts[3]=0;
      //Set new barcode and umi
      curr_barcode = barcode;
    }
    //Add the count to the lowest level counter
    switch(c){
      case 1:
        cnts[0]++;
        break;
      case 2:
        cnts[1]++;
        break;
      case 4:
        cnts[2]++;
        break;
      case 8:
        cnts[3]++;
        break;
      default:
        break;
    }
    p++;
  }
  //Now finalise the last read
  if(n_plp>0){
    //Get the consensus read
    max_obs_reads = -1;
    j = -1;
    for(k=0;k<4;k++){
      if(cnts[k] == max_obs_reads){
        j=-1;
      }
      if(cnts[k] > max_obs_reads){
        max_obs_reads = cnts[k];
        j=k;
      }
    }
    //Add it to the cell level counter
    if(j<0){
       //printf("No consensus allele: %d,%d,%d,%d\n",cnts[0],cnts[1],cnts[2],cnts[3]);
    }else{
      cellCnts[j]++;
    }
    //Print the result (if it's worth printing)
    if(cellCnts[0]+cellCnts[1]+cellCnts[2]+cellCnts[3]>0)
      print_10x_section(output,stats->chr,stats->pos,cellCnts[0],cellCnts[1],cellCnts[2],cellCnts[3],cellCnts[0]+cellCnts[1]+cellCnts[2]+cellCnts[3],curr_barcode);
  }
  return;
}

int bam_access_get_multi_position_base_counts(loci_stats **stats, int stats_count,int is_10x,FILE* output){
	char *region = NULL;
	hts_itr_t *iter = NULL;
	bam1_t* b = NULL;
	bam_plp_t buf;

	//Find start and stop for each contig and retrieve a contig at once
	int start = 0;
	int stop = 0;
	char* this_chr;
	int stop_idx = 0;
	int start_idx = 0;
	while(start_idx<stats_count){
		int i=start_idx;
		stop = stats[i]->pos;
		stop_idx = i;
		this_chr = stats[start_idx]->chr;
		start = stats[start_idx]->pos;
		if(i+1<stats_count){
			i++;
			//Calculate stop of contig
			while(strcmp(this_chr,stats[i]->chr)==0){
				stop = stats[i]->pos;
				stop_idx = i;
				i++;
				if(i==stats_count) break;
			}
		}
		region = malloc((sizeof(char *) * (strlen(this_chr)+1))+sizeof(":")+sizeof("-")+(sizeof(char)*((no_of_digits(start)+no_of_digits(stop))+1)));
		check_mem(region);
		sprintf(region,"%s:%d-%d",this_chr,start,stop);
		// initialize pileup
		buf = bam_plp_init(pileup_func, (void *)fholder);
		bam_plp_set_maxcnt(buf,maxitercnt);
		b = bam_init1();
	  iter = sam_itr_querys(fholder->idx, fholder->head, region);
		int j=start_idx;
		int result;
		const bam_pileup1_t *pl;
		int tid, pos, n_plp = -1;
	  while ((result = sam_itr_next(fholder->in, iter, b)) >= 0) {
       uint8_t *aux_val_bcode;
       uint8_t *aux_val_umi;
       //printf("Got another read \n");
	    if(b->core.qual < min_map_qual || (b->core.flag & exc_flag) || (b->core.flag & inc_flag) != inc_flag) continue;
        //Additional check for properly paired reads - they must be in correct paired end orientation
        if(inc_flag & BAM_FPROPER_PAIR){
          if ((!(b->core.flag & BAM_FMREVERSE) == !(b->core.flag & BAM_FREVERSE))) continue;
        }
       //Extract 10x checks
       if(is_10x){
         aux_val_bcode = bam_aux_get(b,"CB");
         aux_val_umi = bam_aux_get(b,"UB");
         if(!aux_val_bcode || !aux_val_umi)
           continue;
       }
       //printf("Which passed quality checks.\n");
	    bam_plp_push(buf, b);
       //printf("And we pushed it to the buffer.\n");
			while ((pl=bam_plp_next(buf, &tid, &pos, &n_plp)) > 0) {
            //printf("Processing pileup at %d stats at %d\n",pos,stats[j]->pos);
				if(j==stats_count || pos+1>stats[stop_idx]->pos) break;
				while(pos+1>stats[j]->pos){
					if(j==stop_idx) break;
					j++;//WE've finished this position, move on (no cvg?)
				}
				if(pos+1==stats[j]->pos){
              //printf("Doing inner pileup for %d reads.\n",n_plp);
              if(is_10x){
					 pileupCounts10x(pl, n_plp, stats[j],output);
              }else{
					 pileupCounts(pl, n_plp, stats[j]);
              }
				}
            //printf("Processing EOL pileup at %d stats at %d\n",pos,stats[j]->pos);
				if(pos+1>=stats[j]->pos && j==stop_idx) break;
	    }
       //printf("Returning to read loading.\n");
	  }//End of iteration through sam_iter
    check(result>=-1, "Error detected (%d) when trying to iterate through region.",result);
		bam_plp_push(buf, 0); // finalize pileup
		while ((pl=bam_plp_next(buf, &tid, &pos, &n_plp)) > 0) {
			if(j==stats_count || pos+1>stats[stop_idx]->pos) break;
			while(pos+1>stats[j]->pos){
				if(j==stop_idx) break;
				j++;//WE've finished this position, move on (no cvg?)
			}
			if(pos+1==stats[j]->pos){
            //printf("Doing final pileup for %d reads.\n",n_plp);
            if(is_10x){
					 pileupCounts10x(pl, n_plp, stats[j],output);
              }else{
					 pileupCounts(pl, n_plp, stats[j]);
              }
			}
			if(pos+1>=stats[j]->pos && j==stop_idx) break;
		}
		bam_plp_destroy(buf);
		free(region);
		bam_destroy1(b);
		start_idx = stop_idx+1;
	}
	return 0;
	error:
	if(iter) sam_itr_destroy(iter);
	if(b) bam_destroy1(b);
	if(region) free(region);
	return 1;


}

int bam_access_get_position_base_counts(char *chr, int posn, loci_stats *stats,int is_10x,FILE *output){
	char *region = NULL;
	hts_itr_t *iter = NULL;
	bam1_t* b = NULL;
	bam_plp_t buf;
	fholder->stats = stats;

	region = malloc((sizeof(char *) * (strlen(chr)+1))+sizeof(":")+sizeof("-")+(sizeof(char)*((no_of_digits(posn)*2)+1)));
	check_mem(region);
	sprintf(region,"%s:%d-%d",chr,posn,posn);
	fholder->beg = posn;
    fholder->end = posn;

  // initialize pileup
	buf = bam_plp_init(pileup_func, (void *)fholder);
	bam_plp_set_maxcnt(buf,maxitercnt);
  /*
  sam_fetch(fholder->in, fholder->idx, ref, fholder->beg, fholder->end, buf, fetch_algo_func);
  */
  //Replace fetch with iterator for htslib compatibility.
  b = bam_init1();
  iter = sam_itr_querys(fholder->idx, fholder->head, region);
  int result;
  uint8_t *aux_val_bcode;
  uint8_t *aux_val_umi;
  //char *barcode;
  //char *umi;
  while ((result = sam_itr_next(fholder->in, iter, b)) >= 0) {
    if(is_10x){
      aux_val_bcode = bam_aux_get(b,"CB");
      aux_val_umi = bam_aux_get(b,"UB");
      if(!aux_val_bcode || !aux_val_umi){
        continue;
        //printf("Failed to get tags \n");
      }
    }
    if(b->core.qual < min_map_qual || (b->core.flag & exc_flag) || (b->core.flag & inc_flag) != inc_flag) continue;
    bam_plp_push(buf, b);
    //barcode = bam_aux2Z(aux_val_bcode);
    //umi = bam_aux2Z(aux_val_umi);
    //printf("Got tag: bc=%s umi=%s\n",barcode,umi);
  }
  check(result>=-1, "Error detected (%d) when trying to iterate through region.",result);
  sam_itr_destroy(iter);
  bam_plp_push(buf, 0);
  int tid, pos, n_plp = -1;
  const bam_pileup1_t *pil;
  while ( (pil=bam_plp_next(buf, &tid, &pos, &n_plp)) > 0) {
    if((pos+1) != posn) continue;
      if(is_10x){
		  pileupCounts10x(pil, n_plp, fholder->stats,output);
      }else{
		  pileupCounts(pil, n_plp, fholder->stats);
      }
  } //End of iteration through pileup
	//bam_plp_push(buf, 0); // finalize pileup
  bam_plp_destroy(buf);
	free(region);
	bam_destroy1(b);
	return 0;

error:
	//if(region) free(region);
	if(fholder->stats){
		if(fholder->stats->base_counts) free(fholder->stats->base_counts);
		free(fholder->stats);
	}
	if(iter) sam_itr_destroy(iter);
	if(b) bam_destroy1(b);
	if(region) free(region);
	return 1;
}

void bam_access_min_base_qual(int qual){
	min_base_qual = qual;
	return;
}

void bam_access_min_map_qual(int qual){
	min_map_qual = qual;
	return;
}

void bam_access_inc_flag(int inc){
  inc_flag = inc;
  return;
}

void bam_access_exc_flag(int exc){
  exc_flag = exc;
  return;
}
