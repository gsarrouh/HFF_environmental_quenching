# Created on Mon Jul 20 06:47:36 2020
#
################## phot_check.py ##################
#
#
### WHAT THIS PROGRAM IS:
### This short program checks which targets in the "master_cat" catalogue have "good" photometry, i.e. have interpolated values for the U,V,J filters. My understanding is that this would arise if insufficient measurements of sufficient quality are obtained in the (8) HST filters labelled "F***W", where the "***" is a numerical designation, e.g. F160W or F814W. it also identifies galaxies which have been flagged using the "use_phot" flag system as described in Shipley et al 2018, S3.9.
#
## It produces a summary table showing the result
#
## modules
import numpy as np
#
use_phot_counter0 = np.array([0]*6)
use_phot_counter1 = np.array([0]*6)

count_bcg = np.array([0]*6)
uvj_counter = np.array([0]*6)
F435W_counter = np.array([0]*6)
F606W_counter = np.array([0]*6)
F814W_counter = np.array([0]*6)
F105W_counter = np.array([0]*6)
F125W_counter = np.array([0]*6)
F140W_counter = np.array([0]*6)
F160W_counter = np.array([0]*6)
neg_99_counter = np.array([[0]*6]*8)
#
for counter in range(len(master_cat)):
    for cluster in range(len(z_cluster)):
        if master_cat['cluster'][counter] == (cluster+1):         # track by cluster
            # track objects w/ bad "use_phot" flag
            if master_cat['use_phot'][counter] == 0:
                use_phot_counter0[cluster]+=1
            elif master_cat['use_phot'][counter] == 1:
                use_phot_counter1[cluster]+=1
            # track objects designated as bcg's that were modelled out
            if master_cat['bandtotal'][counter] == 'bcg':
                count_bcg[cluster]+=1
            # track objects w/ bad EAZY output photometric rest-frame estimates
            if master_cat['u'][counter] < 0 and master_cat['v'][counter] < 0 and master_cat['j'][counter] < 0:
                uvj_counter[cluster]+=1
                if master_cat['u'][counter] == -99 and master_cat['v'][counter] == -99 and master_cat['j'][counter] == -99:
                    neg_99_counter[0][cluster]+=1
            if master_cat['f_F435W'][counter] < 0:
                F435W_counter[cluster]+=1
                if master_cat['f_F435W'][counter] ==-99:
                    neg_99_counter[1][cluster]+=1
            # bad measurement in F606W filter
            if master_cat['f_F606W'][counter] < 0:
                F606W_counter[cluster]+=1
                if master_cat['f_F606W'][counter] ==-99:
                    neg_99_counter[2][cluster]+=1
            # bad measurement in F814W filter
            if master_cat['f_F814W'][counter] < 0:
                F814W_counter[cluster]+=1
                if master_cat['f_F814W'][counter] ==-99:
                    neg_99_counter[3][cluster]+=1
            if master_cat['f_F105W'][counter] < 0:
                F105W_counter[cluster]+=1
                if master_cat['f_F105W'][counter] ==-99:
                    neg_99_counter[4][cluster]+=1
            # bad measurement in F135W filter
            if master_cat['f_F125W'][counter] < 0:
                F125W_counter[cluster]+=1
                if master_cat['f_F125W'][counter] ==-99:
                    neg_99_counter[5][cluster]+=1
            # bad measurement in F140W filter
            if master_cat['f_F140W'][counter] < 0:
                F140W_counter[cluster]+=1
                if master_cat['f_F140W'][counter] ==-99:
                    neg_99_counter[6][cluster]+=1
            # bad measurement in F160W filter
            if master_cat['f_F160W'][counter] < 0:
                F160W_counter[cluster]+=1
                if master_cat['f_F160W'][counter] ==-99:
                    neg_99_counter[7][cluster]+=1
                
#
## Summarize good_phot stats in table
good_phot_names = Column(['Total','flux(U,V,&J) < 0','flux(ACS_F435W) < 0','flux(ACS_F606W) < 0','flux(ACS_F814W) < 0','flux(WFC3_F105W) < 0','flux(WFC3_F125W) < 0','flux(WFC3_F140W) < 0','flux(WFC3_F160W) < 0','bCGs modelled out'],name='Property')
col_names = cluster_names
good_phot0 = Column([np.sum([both,phot_only,spec_only,no_data,stars_sub]),np.sum(uvj_counter),np.sum(F435W_counter),np.sum(F606W_counter),np.sum(F814W_counter),np.sum(F105W_counter),np.sum(F125W_counter),np.sum(F140W_counter),np.sum(F160W_counter),np.sum(count_bcg)],name='Total')  # total column
good_phot_stats = Table([good_phot_names,good_phot0])
for ii in range(len(num_BCG)):
    good_phot_col = Column([np.sum([both[ii],phot_only[ii],spec_only[ii],no_data[ii],stars_sub[ii]]),uvj_counter[ii],F435W_counter[ii],F606W_counter[ii],F814W_counter[ii],F105W_counter[ii],F125W_counter[ii],F140W_counter[ii],F160W_counter[ii],count_bcg[ii]],name=col_names[ii])               # cluster columns
    good_phot_stats.add_column(good_phot_col)
#
#
## Summarize neg99_phot stats in table
neg99_phot_names = Column(['Total','flux(U,V,&J) = -99','flux(ACS_F435W) = -99','flux(ACS_F606W) = -99','flux(ACS_F814W) = -99','flux(WFC3_F105W) = -99','flux(WFC3_F125W) = -99','flux(WFC3_F140W) = -99','flux(WFC3_F160W) = -99','bCGs modelled out'],name='Property')
col_names = cluster_names
neg99_phot0 = Column([np.sum([both,phot_only,spec_only,no_data,stars_sub]),np.sum(neg_99_counter[0]),np.sum(neg_99_counter[1]),np.sum(neg_99_counter[2]),np.sum(neg_99_counter[3]),np.sum(neg_99_counter[4]),np.sum(neg_99_counter[5]),np.sum(neg_99_counter[6]),np.sum(neg_99_counter[7]),np.sum(count_bcg)],name='Total')  # total column
neg99_phot_stats = Table([neg99_phot_names,neg99_phot0])
for ii in range(len(num_BCG)):
    neg99_phot_col = Column([np.sum([both[ii],phot_only[ii],spec_only[ii],no_data[ii],stars_sub[ii]]),neg_99_counter[0][ii],neg_99_counter[1][ii],neg_99_counter[2][ii],neg_99_counter[3][ii],neg_99_counter[4][ii],neg_99_counter[5][ii],neg_99_counter[6][ii],neg_99_counter[7][ii],count_bcg[ii]],name=col_names[ii])               # cluster columns
    neg99_phot_stats.add_column(neg99_phot_col)
#
#
## Summarize difference stats in table
diff_phot_names = Column(['flux(U,V,&J)','flux(ACS_F435W)','flux(ACS_F606W)','flux(ACS_F814W)','flux(WFC3_F105W)','flux(WFC3_F125W)','flux(WFC3_F140W)','flux(WFC3_F160W)'],name='Property')
col_names = cluster_names
diff_phot0 = Column([np.sum(uvj_counter)-np.sum(neg_99_counter[0]),np.sum(F435W_counter)-np.sum(neg_99_counter[1]),np.sum(F606W_counter)-np.sum(neg_99_counter[2]),np.sum(F814W_counter)-np.sum(neg_99_counter[3]),np.sum(F105W_counter)-np.sum(neg_99_counter[4]),np.sum(F125W_counter)-np.sum(neg_99_counter[5]),np.sum(F140W_counter)-np.sum(neg_99_counter[6]),np.sum(F160W_counter)-np.sum(neg_99_counter[7])],name='Total')  # total column
diff_phot_stats = Table([diff_phot_names,diff_phot0])
for ii in range(len(num_BCG)):
    diff_phot_col = Column([uvj_counter[ii]-neg_99_counter[0][ii],F435W_counter[ii]-neg_99_counter[1][ii],F606W_counter[ii]-neg_99_counter[2][ii],F814W_counter[ii]-neg_99_counter[3][ii],F105W_counter[ii]-neg_99_counter[4][ii],F125W_counter[ii]-neg_99_counter[5][ii],F140W_counter[ii]-neg_99_counter[6][ii],F160W_counter[ii]-neg_99_counter[7][ii]],name=col_names[ii])               # cluster columns
    diff_phot_stats.add_column(diff_phot_col)
#
#
#
print('\nSummary Table: General photometry (FLUX < 0) stats\n%s'%good_phot_stats)
print('\nSummary Table: General photometry (FLUX == -99) stats\n%s'%neg99_phot_stats)
print('\nDifference Table: (flux < 0) - (flux == -99) \n%s'%diff_phot_stats)
print('\nNOTE: # of bCGs identified as having been modelled out: %s'%np.sum(count_bcg))
#
#
#
#
## Summarize neg99_phot stats in table
use_phot_names = Column(['use_phot == 0','use_phot == 1'],name='Property')
col_names = cluster_names
use_phot0 = Column([np.sum(use_phot_counter0),np.sum(use_phot_counter1)],name='Total')  # total column
use_phot_stats = Table([use_phot_names,use_phot0])
for ii in range(len(num_BCG)):
    use_phot_col = Column([use_phot_counter0[ii],use_phot_counter1[ii]],name=col_names[ii])               # cluster columns
    use_phot_stats.add_column(use_phot_col)
#
print('\nSummary Table: "use_phot" stats\n%s'%use_phot_stats)


#
#
#
###################     PROGRAM END     ###################