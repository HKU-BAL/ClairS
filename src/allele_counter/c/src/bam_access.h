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

#ifndef _bam_access_h
#define _bam_access_h

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <htslib/sam.h>

typedef struct loci_stats{
	int *base_counts;
	char *chr;
	int pos;
	char allele_A;
	char allele_B;
} loci_stats;

typedef struct file_holder{
	int beg, end;
	htsFile *in;
	hts_idx_t *idx;
	loci_stats *stats;
	bam_hdr_t *head;
} file_holder;

void bam_access_min_base_qual(int qual);

void bam_access_min_map_qual(int qual);

void bam_access_inc_flag(int inc);

void bam_access_exc_flag(int exc);

int bam_access_openhts(char *hts_file, char *ref_file);

int bam_access_get_position_base_counts(char *chr, int pos, loci_stats *stats,int is_10x,FILE *output);

int bam_access_get_multi_position_base_counts(loci_stats **stats, int stats_count,int is_10x,FILE *output);

void bam_access_closehts();

int readCompare(const void *r1,const void *r2);

#endif
