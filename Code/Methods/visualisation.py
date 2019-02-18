# DRG Pom visualisation module, adapted from work by Alonso and Marder, eLife, 2019.



from pylab import *

def plot_currentscape(voltage, currents, figsize=(3,4)):	
	# make a copy of currents
	# CURRENTSCAPE CALCULATION STARTS HERE. 
	curr=array(currents)	
	cpos= curr.copy()
	cpos[curr<0]=0
	cneg= curr.copy()
	cneg[curr>0]=0

	normapos = sum(abs(array(cpos)),axis=0)
	normaneg = sum(abs(array(cneg)),axis=0)
	npPD=normapos
	nnPD=normaneg
	cnorm=curr.copy()
	cnorm[curr>0]=(abs(curr)/normapos)[curr>0]
	cnorm[curr<0]=-(abs(curr)/normaneg)[curr<0]

	resy=1000
	impos=zeros((resy,shape(cnorm)[-1])) 
	imneg=zeros((resy,shape(cnorm)[-1])) 

	times=arange(0,shape(cnorm)[-1])
	for t in times:
	    lastpercent=0
	    for numcurr, curr in enumerate(cnorm):
	        if(curr[t]>0):
	            percent = int(curr[t]*(resy))   
	            impos[lastpercent:lastpercent+percent,t]=numcurr
	            lastpercent=lastpercent+percent        
	for t in times:
	    lastpercent=0
	    for numcurr, curr in enumerate(cnorm):
	        if(curr[t]<0):
	            percent = int(abs(curr[t])*(resy))   
	            imneg[lastpercent:lastpercent+percent,t]=numcurr
	            lastpercent=lastpercent+percent        
	im0= vstack((impos,imneg))
	# CURRENTSCAPE CALCULATION ENDS HERE. 

	#PLOT CURRENTSCAPE
	fig = figure(figsize=figsize)

	#PLOT VOLTAGE TRACE
	xmax=len(voltage)
	swthres=-50        
	ax=subplot2grid((7,1),(0,0),rowspan=2)	
	t=arange(0,len(voltage))
	plot(t, voltage, color='black',lw=1.)
	plot(t,ones(len(t))*swthres,ls='dashed',color='black',lw=0.75)
	vlines(1,-50,-20,lw=1)
	ylim(-75,70)
	xlim(0,xmax)
	axis('off')         

	#PLOT TOTAL INWARD CURRENT IN LOG SCALE
	ax=subplot2grid((7,1),(2,0),rowspan=1)
	fill_between(arange(len((npPD))),(npPD),color='black')
	plot(5.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(50.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(500.*ones(len(nnPD)),color='black', ls=':',lw=1)
	yscale('log')
	ylim(0.01,1500)
	xlim(0,xmax)
	axis('off') 

	#PLOT CURRENT SHARES
	elcolormap='Set1'
	ax=subplot2grid((7,1),(3,0),rowspan=3)
	imshow(im0[::1,::1],interpolation='nearest',aspect='auto',cmap=elcolormap)
	ylim(2*resy,0)
	plot(resy*ones(len(npPD)),color='black',lw=2)
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	plt.gca().yaxis.set_major_locator(plt.NullLocator())
	xlim(0,xmax)
	clim(0,8)
	axis('off') 

	#PLOT TOTAL OUTWARD CURRENT IN LOG SCALE
	ax=subplot2grid((7,1),(6,0),rowspan=1)
	fill_between(arange(len((nnPD))),(nnPD),color='black')
	plot(5.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(50.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(500.*ones(len(nnPD)),color='black', ls=':',lw=1)
	yscale('log')
	ylim(1500,0.01)
	xlim(0,xmax)
	axis('off') 
	subplots_adjust(wspace=0, hspace=0)
	return fig
plotCurrentscape = plot_currentscape

def plot_voltage_distributions(Vdist): 
    im=Vdist
    fig = figure()
    cmmap='Greys'		
    #proper visualization requires that we take the logarithm (we add +1 to avoid NaN)       
    imshow((log10(im+1)),aspect='auto', cmap=cmmap, extent=[1,0,-75,30], interpolation='nearest')
    xlim(1,0)
    ylim(-75,30)
    clim(1,5)
    axis('off')     
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)    
    return fig
plotVoltage_distributions = plot_voltage_distributions

def plot_voltage_distributions_enhance_edges(Vdist): 
    im=Vdist
    fig = figure(figsize=(3,6))    
    #choose a colormap (see also gnuplot1, gnuplot2, helix)
    cmmap='hot'
    #proper visualization requires that we take the logarithm     
    a=log10(im+1)
    #We smooth the distribution using a convolution with the identity. This does not have a big effect. 
    filt=ones(3)/3.
    r=np.apply_along_axis(lambda m: np.convolve(m, filt, mode='same'), axis=0, arr=a)
    #We enhance edges by performing a derivative (along the V-axis). 
    imm=abs(diff(r,axis=0))
    imshow(imm,aspect='auto', cmap=cmmap, extent=[1,0,-75,30], interpolation='nearest')     
    xlim(1,0)
    ylim(-75,30)
    axis('off')     
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)    
    clim(0.0,0.3)
    return fig
plotVoltageDistributionsEnhanceEdges = plot_voltage_distributions_enhance_edges

def plot_current_shares_distributions(current_share_dist): 
    im=current_share_dist
    fig = figure()
    cmmap='gnuplot2'		
    #proper visualization requires that we take the logarithm (we add +1 to avoid NaN)       
    imshow(log10(im+1),aspect='auto', cmap=cmmap, extent=[1,0,0,1], interpolation='nearest')
    percents=linspace(1,0,101)
    plot(percents,ones(len(percents))*(0.25), color = 'white',ls=':')    
    plot(percents,ones(len(percents))*(0.5), color = 'white',ls=':')    
    plot(percents,ones(len(percents))*(0.75), color = 'white',ls=':')   
    xlim(1,0)
    ylim(0,1)
    clim(1,5)
    axis('off')     
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)    
    return fig    
plotCurrentSharesDistributions = plot_current_shares_distributions

def plot_pie_chart(currents, title='', legend=['Nav17','Nav18','Nav19','Kdr','KA','KM','HCN'], autopct=None):
    """
    Plot pie chart
    
    Inputs:
    currents - an array of currents, each array element should be a float - not a vector of currents like the original pie chart function in Alonso Marder
    title - title string for plot
    legend - ordered list of current names
    """
    # Plot pie charts and batons
    num_currents = currents.shape[0]+1
    color_set = plt.cm.Set1(np.arange(num_currents)/num_currents)

    # Positive currents
    c_pos= currents.copy()
    c_pos[currents<0]=0

    # Negative currents:w

    c_neg=currents.copy()
    c_neg[currents>0]=0

    # Sum of all positive currents and all negative currents
    norm_pos = np.sum(abs(np.array(c_pos)),axis=0)
    norm_neg = np.sum(abs(np.array(c_neg)),axis=0)

    # Get normalised currents for positive and negative currents
    c_norm = currents.copy()
    c_norm[currents>0] = (abs(currents)/norm_pos)[currents>0]
    c_norm[currents<0] = -(abs(currents)/norm_neg)[currents<0]   

    c_norm_pos = np.zeros(np.shape(c_norm))
    c_norm_pos[c_norm>0] = c_norm[c_norm>0]

    c_norm_neg = np.zeros(np.shape(c_norm))
    c_norm_neg[c_norm<0] = c_norm[c_norm<0]
   
    # n not needed as function takes one set of currents, not vectors of currents
    sizes_pos = c_norm_pos
    sizes_neg = abs(c_norm_neg)
    #sizes=sizes_pos/2
    #print(sum(sizes))

    fig1, ax1 = plt.subplots()    
    ax1.pie(sizes_neg[::1]/2, colors=color_set[::1], autopct=autopct, shadow=False, startangle=0, counterclock=False)
    ax1.axis('equal') 


    ax1.pie(sizes_pos[::1]/2,colors=color_set[::1], autopct=autopct, shadow=False, startangle=180,counterclock=False)
    ax1.axis('equal') 

    plt.legend(legend)
    plt.plot(np.linspace(-1.1,1.1,10),np.zeros(10), lw=5,color='black')
    plt.xlim(-1.05,1.05)
    plt.ylim(-1.05,1.05)
    #plt.title('{:.3g} ms'.format(n*delta_t))
    plt.title(title)

    plt.axis('off')         

" ---- Utilities ---- "

def format_trace_for_currentscape(trace, current_names=None, return_current_names=False):
    
    voltage_trace = trace['v']
    if current_names == None:
        current_names = [key for key in trace.keys() if key not in ['t', 'v']]
        
    currents = np.zeros( (len(current_names), len(voltage_trace)) ) # Array has one row per current and same length as voltage

    for i, current_name in enumerate(current_names):
        # Combine the currents if there's more than one
        for current_component in trace[current_name]:
            currents[i] = currents[i] + trace[current_name][current_component]
    
    if return_current_names:
        return voltage_trace, currents, current_names
    else:
        return voltage_trace, currents


	

