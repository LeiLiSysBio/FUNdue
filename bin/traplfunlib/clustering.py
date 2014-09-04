import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pylab
import matplotlib as mpl
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy as np
import string
import time
import sys, os
import getopt

################# Perform the hierarchical clustering #################
class cluster_analysis(object):   
    def __init__(self,filename, expressionpath, rowmethod, rowmatrix, colmethod, colmatrix, color, idcol):
        self._filename = filename
        self._expressionpath = expressionpath
        self._rowmethod = rowmethod
        self._rowmatrix = rowmatrix
        self._colmethod = colmethod
        self._colmatrix = colmatrix
        self._color = color
        self._idcol = idcol
                   
    def _importfromData(self,filename):
        matrix=[]
        row_header=[]
        first_row=True
    
        if '/' in filename:
            dataset_name = filename.split('/')[-1][:-4]
        else:
            dataset_name = filename.split('\\')[-1][:-4]
    
        for line in open(filename,'r'):
            t = line[:-1].split('\t') ### remove end-of-line character - file is tab-delimited
            if first_row:
                column_header = t[1:]
                first_row=False
            else:
                if ' ' not in t and '' not in t: ### Occurs for rows with missing data
                    s = list(map(float,t[1:]))
                    if (abs(max(s)-min(s)))>0:
                        matrix.append(s)
                        row_header.append(t[0])   
                        
        return np.array(matrix), column_header, row_header
    
    def _heatmap(self, filematrix, row_header, column_header, row_method,
            column_method, row_metric, column_metric,
            color_gradient, filename):
        
        """x is an m by n ndarray, m observations, n genes"""
        ### Define the color gradient to use based on the provided name
        n = len(filematrix[0]) 
        m = len(filematrix)
        if color_gradient == 'red_white_blue':
            cmap=pylab.cm.bwr
        if color_gradient == 'yellow_red':
            cmap=pylab.cm.YlOrRd
        if color_gradient == 'red_black_green':
            cmap=self._RedBlackGreen()
        if color_gradient == 'seismic':
            cmap=pylab.cm.seismic
        if color_gradient == 'green_white_purple':
            cmap=pylab.cm.PiYG_r
        if color_gradient == 'coolwarm':
            cmap=pylab.cm.coolwarm
        ### Scale the max and min colors so that 0 is white/black
        vmin=filematrix.min()
        vmax=filematrix.max()
        vmax = max([vmax,abs(vmin)])
        vmin = vmax*-1
        norm = mpl.colors.Normalize(vmin/2, vmax/2) ### adjust the max and min to scale these colors
        ### Scale the Matplotlib window size
        default_window_hight = 8.5
        default_window_width = 12
        fig = pylab.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
        color_bar_w = 0.015 ### Sufficient size to show
    
        ## calculate positions for all elements
        # ax1, placement of dendrogram 1, on the left of the heatmap
        [ax1_x, ax1_y, ax1_w, ax1_h] = [0.18,0.22,0.1,0.6]   ### The second value controls the position of the matrix relative to the bottom of the view
        width_between_ax1_axr = 0.004
        height_between_ax1_axc = 0.004 ### distance between the top color bar axis and the matrix
    
        # axr, placement of row side colorbar
        [axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w,0.6] ### second to last controls the width of the side color bar - 0.015 when showing
        axr_x = ax1_x + ax1_w + width_between_ax1_axr
        axr_y = ax1_y; axr_h = ax1_h
        width_between_axr_axm = 0.004
    
        # axc, placement of column side colorbar
        [axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.5,color_bar_w] ### last one controls the hight of the top color bar - 0.015 when showing
        axc_x = axr_x + axr_w + width_between_axr_axm
        axc_y = ax1_y + ax1_h + height_between_ax1_axc
        height_between_axc_ax2 = 0.004
    
        # axm, placement of heatmap for the data matrix
        [axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5]
        axm_x = axr_x + axr_w + width_between_axr_axm
        axm_y = ax1_y; axm_h = ax1_h
        axm_w = axc_w
    
        # ax2, placement of dendrogram 2, on the top of the heatmap
        [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.08] ### last one controls hight of the dendrogram
        ax2_x = axr_x + axr_w + width_between_axr_axm
        ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
        ax2_w = axc_w
    
        # axcb - placement of the color legend
        [axcb_x, axcb_y, axcb_w, axcb_h] = [0.07,0.88,0.18,0.09]
    
        # Compute and plot top dendrogram
        if column_method.lower() != 'none':
            d2 = dist.pdist(filematrix.T)
            D2 = dist.squareform(d2)
            ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=True)
            Y2 = sch.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
            sch.set_link_color_palette(['black'])
            Z2 = sch.dendrogram(Y2, color_threshold=np.inf)
            ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
            ax2.set_xticks([]) ### Hides ticks
            ax2.set_yticks([])
        else:
            ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data
    
        # Compute and plot left dendrogram.
        if row_method.lower() != 'none':
            d1 = dist.pdist(filematrix)
            D1 = dist.squareform(d1)  # full matrix
            ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True) # frame_on may be False
            Y1 = sch.linkage(D1, method=row_method, metric=row_metric)
            sch.set_link_color_palette(['black'])
            Z1 = sch.dendrogram(Y1, color_threshold=np.inf, orientation='right')
            ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
            ax1.set_xticks([]) ### Hides ticks
            ax1.set_yticks([])
        else:
            ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data
    
        # Plot distance matrix.
        axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])  # axes for the data matrix
        xt = filematrix
        if column_method.lower() != 'none':
            idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
            xt = xt[:,idx2]
            ind2 = [ind2[i] for i in idx2]
        if row_method.lower() != 'none':
            idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
            xt = xt[idx1,:]   # xt is transformed x
            ind1 = [ind1[i] for i in idx1] ### reorder the flat cluster to match the order of the leaves the dendrogram
        im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black
        axm.set_xticks([]) ### Hides x-ticks
        axm.set_yticks([])
    
        # Add text
        new_row_header=[]
        new_column_header=[]
        for i in range(filematrix.shape[0]):
            if row_method.lower() != 'none':
                if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
                    axm.text(filematrix.shape[1]-0.5, i, '  '+row_header[idx1[i]])
                new_row_header.append(row_header[idx1[i]])
            else:
                if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
                    axm.text(filematrix.shape[1]-0.5, i, '  '+row_header[i]) ### When not clustering rows
                new_row_header.append(row_header[i])
        for i in range(filematrix.shape[1]):
            if column_method.lower() != 'none':
                axm.text(i, -0.9, ' '+column_header[idx2[i]], rotation=270, verticalalignment="top") # rotation could also be degrees
                new_column_header.append(column_header[idx2[i]])
            else: ### When not clustering columns
                axm.text(i, -0.9, ' '+column_header[i], rotation=270, verticalalignment="top")
                new_column_header.append(column_header[i])
    
        # Plot colside colors
        # axc --> axes for column side colorbar
        if column_method.lower() != 'none':
            axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
            cmap_c = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
            dc = np.array(ind2, dtype=int)
            dc.shape = (1,len(ind2))
            im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
            axc.set_xticks([]) ### Hides ticks
            axc.set_yticks([])
    
        # Plot rowside colors
        # axr --> axes for row side colorbar
        if row_method.lower() != 'none':
            axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
            dr = np.array(ind1, dtype=int)
            dr.shape = (len(ind1),1)
            #print ind1, len(ind1)
            cmap_r = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
            im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_c)
            axr.set_xticks([]) ### Hides ticks
            axr.set_yticks([])
    
        # Plot color legend
        axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
        cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
    
        if '/' in filename:
            dataset_name = filename.split('/')[-1][:-4]
            root_dir = str.join('/', filename.split('/')[:-1])+'/'
        #else:
            #dataset_name = filename.split('\\')[-1][:-4]
            #root_dir = string.join(string.split(filename,'\\')[:-1],'\\')+'\\'
        filename = self._expressionpath+ 'hierarchical-Clustering-%s-%s_%s.pdf' % (dataset_name,column_metric,row_metric)
        cb.set_label("Enrichment factor")
        #self._exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2)
    
        ### Render the graphic
        if len(row_header)>50 or len(column_header)>50:
            pylab.rcParams['font.size'] = 5
        else:
            pylab.rcParams['font.size'] = 8
    
        pylab.savefig(filename)
        print('The clustering result can be found in',filename)
        filename = filename[:-3]+'png'
        pylab.savefig(filename, dpi=100) #,dpi=200
        pylab.show()
    
    def _getColorRange(x):
        """ Determines the range of colors, centered at zero, for normalizing cmap """
        vmax=x.max()
        vmin=x.min()
        if vmax<0 and vmin<0: direction = 'negative'
        elif vmax>0 and vmin>0: direction = 'positive'
        else: direction = 'both'
        if direction == 'both':
            vmax = max([vmax,abs(vmin)])
            vmin = -1*vmax
            return vmax,vmin
        else:
            return vmax,vmin
    
    def _RedBlackGreen(self):
        cdict = {'red':   ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.1),
                           (1.0, 1.0, 1.0)),
    
                 'blue': ((0.0, 0.0, 0.0),
                           (1.0, 0.0, 0.0)),
    
                 'green':  ((0.0, 0.0, 1.0),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0))
                }
    
        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap
        
    def _exportFlatClusterData(self,filename, new_row_header,new_column_header,xt,ind1,ind2):
        """ Export the clustered results as a text file, only indicating the flat-clusters rather than the tree """
    
        filename = filename.replace('.pdf','.txt')
        export_text = open(filename,'w')
        column_header = str.join('\t', ['UID','row_clusters-flat']+new_column_header)+'\n' ### format column-names for export
        export_text.write(column_header)
        column_clusters = str.join('\t', ['column_clusters-flat','']+ map(str, ind2))+'\n' ### format column-flat-clusters for export
        export_text.write(column_clusters)
    
        ### The clusters, dendrogram and flat clusters are drawn bottom-up, so we need to reverse the order to match
        new_row_header = new_row_header[::-1]
        xt = xt[::-1]
    
        ### Export each row in the clustered data matrix xt
        i=0
        for row in xt:
            export_text.write(str.join('\t',[new_row_header[i],str(ind1[i])]+map(str, row))+'\n')
            i+=1
        export_text.close()
    
        ### Export as CDT file
        filename = string.replace(filename,'.txt','.cdt')
        export_cdt = open(filename,'w')
        column_header = str.join('\t', ['UNIQID','NAME','GWEIGHT']+new_column_header)+'\n' ### format column-names for export
        export_cdt.write(column_header)
        eweight = str.join('\t', ['EWEIGHT','','']+ ['1']*len(new_column_header))+'\n' ### format column-flat-clusters for export
        export_cdt.write(eweight)
    
        ### Export each row in the clustered data matrix xt
        i=0
        for row in xt:
            export_cdt.write(string.join([new_row_header[i]]*2+['1']+map(str, row),'\t')+'\n')
            i+=1
        export_cdt.close()
    
           
    def cluster(self):
        start_time = time.time()
        #euclidean,minkowski,cityblock, sequclidean, cosine, correlaiton, 
        row_method = self._rowmethod
        column_method = self._colmethod
        row_metric = self._rowmatrix
        column_metric = self._colmatrix
        color_gradient = self._color
        matrix, column_header, row_header = self._importfromData(self._filename)
        time_diff = str(round(time.time()-start_time,1))
        if len(matrix)>0:
            try:
                self._heatmap(matrix, row_header, \
                column_header, row_method, \
                column_method, row_metric, \
                column_metric, color_gradient, self._filename)
            except IOError:
                print('Error with clustering encountered')
        
    

    