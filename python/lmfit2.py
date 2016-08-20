
fitacf_types = {
 'radar.revision.major' : 'c',
 'radar.revision.minor' : 'c',
 'origin.code' : 'c',
 'origin.time' : 's',
 'origin.command' : 's',
 'cp' : 'h',
 'stid' : 'h',
 'time.yr' : 'h',
 'time.mo' : 'h',
 'time.dy' : 'h',
 'time.hr' : 'h',
 'time.mt' : 'h',
 'time.sc' : 'h',
 'time.us' : 'i',
 'txpow' : 'h',
 'nave' : 'h',
 'atten' : 'h',
 'lagfr' : 'h',
 'smsep' : 'h',
 'ercod' : 'h',
 'stat.agc' : 'h',
 'stat.lopwr' : 'h',
 'noise.search' : 'f',
 'noise.mean' : 'f',
 'channel' : 'h',
 'bmnum' : 'h',
 'bmazm' : 'f',
 'scan' : 'h',
 'offset' : 'h',
 'rxrise' : 'h',
 'intt.sc' : 'h',
 'intt.us' : 'i',
 'txpl' : 'h',
 'mpinc' : 'h',
 'mppul' : 'h',
 'mplgs' : 'h',
 'nrang' : 'h',
 'frang' : 'h',
 'rsep' : 'h',
 'xcf' : 'h',
 'tfreq' : 'h',
 'tfreq2' : 'h',
 'mxpwr' : 'i',
 'lvmax' : 'i',
 'fitacf.revision.major' : 'i',
 'fitacf.revision.minor' : 'i',
 'combf' : 's',
 'noise.sky' : 'f',
 'noise.lag0' : 'f',
 'noise.vel' : 'f',
 'ptab' : 'h',
 'ltab' : 'h',
 'pwr0' : 'f',
 'slist' : 'h',
 'nlag' : 'h',
 'qflg' : 'c',
 'gflg' : 'c',
 'p_l' : 'f',
 'p_l_e' : 'f',
 'p_s' : 'f',
 'p_s_e' : 'f',
 'v' : 'f',
 'v_e' : 'f',
 'w_l' : 'f',
 'w_l_e' : 'f',
 'w_s' : 'f',
 'w_s_e' : 'f',
}

fit_record_keys = [
	 'radar.revision.major',
	 'radar.revision.minor',
	 'origin.code',
	 'origin.time',
	 'origin.command',
	 'cp',
	 'stid',
	 'time.yr',
	 'time.mo',
	 'time.dy',
	 'time.hr',
	 'time.mt',
	 'time.sc',
	 'time.us',
	 'txpow',
	 'nave',
	 'atten',
	 'lagfr',
	 'smsep',
	 'ercod',
	 'stat.agc',
	 'stat.lopwr',
	 'noise.search',
	 'noise.mean',
	 'channel',
	 'bmnum',
	 'bmazm',
	 'scan',
	 'offset',
	 'rxrise',
	 'intt.sc',
	 'intt.us',
	 'txpl',
	 'mpinc',
	 'mppul',
	 'mplgs',
	 'nrang',
	 'frang',
	 'rsep',
	 'xcf',
	 'tfreq',
	 'mxpwr',
	 'lvmax',
	 'fitacf.revision.major',
	 'fitacf.revision.minor',
	 'combf',
	 'noise.sky',
	 'noise.lag0',
	 'noise.vel',
	 'ptab',
	 'ltab',
	 'pwr0',
	 'slist',
	 'nlag',
	 'qflg',
	 'gflg',
	 'p_l',
	 'p_l_e',
	 'p_s',
	 'p_s_e',
	 'v',
	 'v_e',
	 'w_l',
	 'w_l_e',
	 'w_s',
	 'w_s_e',
	 'sd_l',
	 'sd_s',
	 'sd_phi',
	 'x_qflg',
	 'x_gflg',
	 'x_p_l',
	 'x_p_l_e',
	 'x_p_s',
	 'x_p_s_e',
	 'x_v',
	 'x_v_e',
	 'x_w_l',
	 'x_w_l_e',
	 'x_w_s',
	 'x_w_s_e',
	 'phi0',
	 'phi0_e',
	 'elv',
	 'elv_low',
	 'elv_high',
	 'x_sd_l',
	 'x_sd_s',
	 'x_sd_phi']


#read rawacf file using pydmap and return a list of 
#dictionaries each containing a beam record
def read_file(filename):
    from dmap import parse_dmap_format_from_file

    return parse_dmap_format_from_file(filename)


#write the fitted file using pydmap from a list of 
#dictionaries each containing a beam record
def write_file(filename, records):
    from dmap import RawDmapWrite
    RawDmapWrite(records,filename,ud_types=fitacf_types)


#Estimate teh noise power from the 10 lowest power
#ranges.
def estimate_noise(power):
    import numpy as np
    p = list(power)
    p.sort()
    noise = np.mean(p[0:20])
    return noise


#First order error bars for magnitude and acf fitting
def first_order_errors(power,noise,clutter,nave):
    import numpy as np
    return (power+noise+clutter)/np.sqrt(nave)


#Errors in real and imaginary components of acfs
def acf_error(power,noise,clutter,nave,rho,rho_c):
    import numpy as np
    return (power+noise+clutter)*np.sqrt((1-rho**2)/2. + rho_c**2)/np.sqrt(nave)


#Magnitude component errorbars
def mean_mag_stats_from_gaus(re,im,re_error,im_error,rho_r,rho_i,K,num_samps=10000):
    import numpy as np
    factor = 2*rho_r*rho_i
    m = np.array([re,im])
    c = np.array([[re_error**2,factor*re_error*im_error],
                  [factor*re_error*im_error,im_error**2]])
    temp = np.random.multivariate_normal(m,c,size=num_samps)
    x = temp[:,0]
    y = temp[:,1]
    r = np.sqrt(x**2 + y**2)

    if r is not None:
        mean = np.mean(r)
        var = np.var(r)
    else:
        mean = np.nan
        var = np.nan
    return mean, var

#Figure out which lags are bad due to transmitter blanking
#the receivers.
def determine_tx_blanked(nrg,ltab,tp,tau,tfr,gate):
    import numpy as np

    #Calculate the lags and the pulse table
    lags = []
    ptab = []
    for pair in ltab:
        lags.append(pair[1]-pair[0])
        ptab.extend(pair)
    lags=list(set(lags))
    ptab=list(set(ptab))
    lags.sort()  
    ptab.sort()

    txs_in_lag={}
    for lag in lags:
        txs_in_lag[lag]=[]

    #number of ranges per tau
    tp_in_tau = tau/tp

    #determine which range gates correspond to blanked samples
    tx_times=[p*tp_in_tau for p in ptab]
    blanked_samples=[]
    for tx in tx_times:
        blanked_samples.extend([tx,tx+1])

    for lag in lags:
        for pair in ltab:
            if (pair[1] - pair[0] == lag):
                #which samples are we using?
                S1=tp_in_tau*pair[0]+gate + 2.*tfr/(tp)
                S2=tp_in_tau*pair[1]+gate + 2.*tfr/(tp)
                br=[]
                #check to see if the samples are blanked or not
                if S1 in blanked_samples:
                    br.append(S1)
                if S2 in blanked_samples:
                    br.append(S2)
                txs_in_lag[lag].extend(br)

    return txs_in_lag

#Estimate upper limit of self clutter using the MPSE of 
#Reimer and Hussey 2015 (radio science)
def estimate_selfclutter(nrg,ltab,tp,tau,tfr,gate,power):
    import numpy as np

    #Calculate the lags and the pulse table
    lags = []
    ptab = []
    for pair in ltab:
        lags.append(pair[1]-pair[0])
        ptab.extend(pair)
    lags=list(set(lags))
    ptab=list(set(ptab))
    lags.sort()  
    ptab.sort()

    ranges_in_lag={}
    samples_in_lag={}
    for lag in lags:
        ranges_in_lag[lag]=[]
        samples_in_lag[lag]=[]

    #number of ranges per tau
    tp_in_tau = tau/tp

    for lag in lags:
        for pair in ltab:
            if (pair[1] - pair[0] == lag):
                #which samples are we using?
                S1=tp_in_tau*pair[0]+gate + 2.*tfr/(tp)
                S2=tp_in_tau*pair[1]+gate + 2.*tfr/(tp)
                samples_in_lag[lag].append([S1, S2])
                #which pulses were transmitted before the samples were recorded
                P1=[pulse for pulse in ptab if pulse*tp_in_tau <= S1]
                P2=[pulse for pulse in ptab if pulse*tp_in_tau <= S2]
                P1.sort()
                P2.sort()
                #what ranges are the previously transmitted pulses coming from
                r1 = [int(S1 - p*tp_in_tau - 2.*tfr/tp) for p in P1]
                r2 = [int(S2 - p*tp_in_tau - 2.*tfr/tp) for p in P2]
                r1.sort()
                r2.sort()
                ranges_in_lag[lag].append([r1, r2])

    cri_power = np.zeros((lags[-1]+1,))

    for l,lag in enumerate(lags):
        cris=[]

        S1=samples_in_lag[lag][0][0]
        S2=samples_in_lag[lag][0][1]
   
        #reimer method where the power profile is used to estimate the self-clutter
        #so that C is the geometric mean of the power from interfering ranges
        for range_pairs in ranges_in_lag[lag]:
            interfering_ranges = []
            #only use ranges that we have pwr0 for
            interfering_S1 = np.array([int(x) for x in range_pairs[0] if ((x != gate) & (x >= 0.) & (x < nrg-1))])
            interfering_S2 = np.array([int(x) for x in range_pairs[1] if ((x != gate) & (x >= 0.) & (x < nrg-1))])

            #First term in the summation for the self-clutter estimate
            term1 = 0
            if interfering_S1.size > 0:
                for ind in interfering_S1:
                    term1 += np.sqrt(power[gate]*power[ind])

            #Second term in the summation for the self-clutter estimate
            term2 = 0
            if interfering_S2.size > 0:
                for ind in interfering_S2:
                    term2 += np.sqrt(power[gate]*power[ind])

            #Third term in the summation for the self-clutter estimate
            term3 = 0
            if (interfering_S1.size > 0 and interfering_S2.size > 0):
                for i in interfering_S1:
                    for ind in interfering_S2:
                        term3 += np.sqrt(power[i]*power[ind])
            cris.append(term1 + term2 + term3)
        cri_power[lag] += cris[:]

    return cri_power


#Build the initial fit_record to write to file
def build_fit_record(record):
    import numpy as np

    fit_record = dict()
    key_list = ['radar.revision.major','radar.revision.minor','origin.code',
                'origin.time','origin.command','cp','stid','time.yr',
                'time.mo','time.dy','time.hr','time.mt','time.sc','time.us',
                'txpow','nave','atten','lagfr','smsep','ercod','stat.agc',
                'stat.lopwr','noise.search','noise.mean','channel','bmnum',
                'bmazm','scan','offset','rxrise','intt.sc','intt.us','txpl',
                'mpinc','mppul','mplgs','nrang','frang','rsep','xcf','tfreq',
                'mxpwr','lvmax','fitacf.revision.major','fitacf.revision.minor',
                'combf','noise.sky','noise.lag0','noise.vel','ptab','ltab','pwr0']
    for key in key_list:
        if key in record.keys():
            fit_record[key] = record[key]

    return fit_record


def acf_residual(params, time, re, im, re_error, im_error, lamda):
    import numpy as np
    pr = params['power']
    wd = params['width']
    vd = params['velocity']

    model1 = pr*np.exp(-time*wd*2*np.pi/lamda)*np.cos(4.0*np.pi*vd*time/lamda)
    model2 = pr*np.exp(-time*wd*2*np.pi/lamda)*np.sin(4.0*np.pi*vd*time/lamda)

    return np.sqrt( ( (re-model1)/re_error )**2 + ( (im-model2)/im_error )**2 )


#With a record dictionary from pydmap fit the contents
def lmfit2(record):
    import numpy as np
    from datetime import datetime
    from lmfit import Minimizer, Parameters

    now = datetime.now()

    #first get globally used values
    tfreq = record['tfreq'] * 1000.0
    mpinc = record['mpinc'] / 1000000.0
    smsep = record['smsep'] / 1000000.0
    lagfr = record['lagfr'] / 1000000.0
    nrang = record['nrang']
    ltab = record['ltab']
    mplgs = record['mplgs']
    nave = record['nave']
    acfd = record['acfd']
    lag0_power = np.array(record['pwr0'])

    #ignore the second lag 0 in the lag table
    ltab = ltab[0:-1]

    fit_record = build_fit_record(record)

    c = 299792458.0
    lamda = c/tfreq
    k = 2*np.pi/lamda
    nyquist_velocity = lamda/(4.*mpinc)

    lags=[]
    for pair in ltab:
        lags.append(pair[1]-pair[0])
    lags.sort()
    max_lag = lags[-1]

    t = np.array(lags)*mpinc
    ranges = np.arange(0,nrang)

    #Estimate the noise
    noise = estimate_noise(lag0_power)

    #Setup the fitted parameter lists
    fitted_power = list()
    fitted_width = list()
    fitted_vels = list()
    fitted_phis = list()
    fitted_power_e = list()
    fitted_width_e = list()
    fitted_vels_e = list()
    fitted_phis_e = list()
    slist = list()

    #next, iterate through all range gates and do fitting
    for gate in ranges:
        print gate
        re = [x[0] for x in acfd[gate]]
        im = [x[1] for x in acfd[gate]]
        time = list(t)

        #estimate the upper limit of the self clutter
        clutter = list(estimate_selfclutter(nrang,ltab,smsep,mpinc,lagfr/2,gate,lag0_power)[lags])

        #find lags blanked due to Tx and identify "good" lags
        blanked = determine_tx_blanked(nrang,ltab,smsep,mpinc,lagfr/2,gate)
        blank_lags = [0 if len(blanked[x])==0 else 1 for x in blanked]
        
        #don't include any lags that are blanked
        j=0
        for i,bl in enumerate(blank_lags):
            if (bl == 1):
                re.pop(i-j)
                im.pop(i-j)
                time.pop(i-j)
                clutter.pop(i-j)
                j += 1
        re = np.array(re)
        im = np.array(im)
        time = np.array(time)
        clutter = np.array(clutter)

        #Now fit each acf using first order errors
        first_error = first_order_errors(lag0_power[gate],noise,clutter,nave)

        #set up the fitter and fit for the first time
        init_vels = np.linspace(-nyquist_velocity/2.,nyquist_velocity/2.,num=30)
        outs = list()

        for vel in init_vels:
            params = Parameters()
            params.add('power', value=lag0_power[gate])
            params.add('width', value=200.0,min=-100) #Need minimum, to stop magnitude model from diverging to infinity
            params.add('velocity', value=vel, min=-nyquist_velocity/2., max=nyquist_velocity/2.)

            minner = Minimizer(acf_residual, params, fcn_args=(time, re, im, first_error, first_error, lamda))
            outs.append(minner.minimize())

        chi2 = np.array([out.chisqr for out in outs])
        ind = np.where(chi2 == np.min(chi2))[0]

        if (ind.size != 1):
            print "SOMETHING WEIRD IS HAPPENING"
        else:
            ind = ind[0]

        pwr_fit = outs[ind].params['power'].value
        wid_fit = outs[ind].params['width'].value
        vel_fit = outs[ind].params['velocity'].value

        #Now get proper errorbars using fitted parameters and model
        acf_model = pwr_fit*np.exp(-time*2.*np.pi*wid_fit/lamda)*np.exp(1j*4.*np.pi*vel_fit*time/lamda)
        mag_model = np.abs(acf_model)
        rho_re = np.cos(4.*np.pi*vel_fit*time/lamda)
        rho_im = np.sin(4.*np.pi*vel_fit*time/lamda)
        rho = mag_model/mag_model[0]
        for i in range(len(rho)):
            if (rho[i] > 0.999):
                rho[i] = 0.999
            rho[i] = rho[i] * pwr_fit / (pwr_fit + noise + clutter[i])

        re_error = acf_error(pwr_fit,noise,clutter,nave,rho,rho_re)
        im_error = acf_error(pwr_fit,noise,clutter,nave,rho,rho_im)

        #Now second LMFIT
        outs2 = list()

        for vel in init_vels:
            params = Parameters()
            params.add('power', value=pwr_fit)
            params.add('width', value=wid_fit,min=-100)
            params.add('velocity', value=vel, min=-nyquist_velocity/2., max=nyquist_velocity/2.)

            minner = Minimizer(acf_residual, params, fcn_args=(time, re, im, re_error, im_error, lamda))
            outs2.append(minner.minimize())

        chi2 = np.array([out.chisqr for out in outs2])
        ind = np.where(chi2 == np.min(chi2))[0]

        if (ind.size != 1):
            print "SOMETHING WEIRD IS HAPPENING"
        else:
            ind = ind[0]

        pwr_fit = outs2[ind].params['power'].value
        wid_fit = outs2[ind].params['width'].value
        vel_fit = outs2[ind].params['velocity'].value
        
        pwr_e = outs2[ind].params['power'].stderr
        wid_e = outs2[ind].params['width'].stderr
        vel_e = outs2[ind].params['velocity'].stderr

        #Now save fitted quantities into array
        slist.append(gate)
        fitted_power.append(pwr_fit)
        fitted_width.append(wid_fit)
        fitted_vels.append(vel_fit)
        fitted_power_e.append(pwr_e)
        fitted_width_e.append(wid_e)
        fitted_vels_e.append(vel_e)

    print "It took "+str((datetime.now()-now).total_seconds())+" to fit one beam."

    #set ground scatter flags
    gflg = list()
    p_l = list()
    p_l_e = list()
    for i in range(len(slist)):
        if (np.abs(fitted_vels[i])-(30.-1./3.*np.abs(fitted_width[i])) < 0.):
            gflg.append(1)
        else:
            gflg.append(0)
        p_l.append(10.0*np.log10(fitted_power[i]/noise))
        p_l_e.append(10.0*np.log10((fitted_power_e[i]+fitted_power[i])/noise)-10.0*np.log10(fitted_power[i]/noise))

    #construct the fitted data dictionary that will be written to the fit file
    fit_record['slist'] = np.array(slist,dtype=np.int16)
    fit_record['nlag'] = mplgs * np.ones(len(slist),dtype=np.int16)
    fit_record['qflg'] = [1]*len(slist)
    fit_record['gflg'] = gflg
    fit_record['p_l'] = np.array(p_l,dtype=np.float32)
    fit_record['p_l_e'] = np.array(p_l_e,dtype=np.float32)
    fit_record['noise.sky'] = noise
    fit_record['noise.search'] = noise
    fit_record['noise.mean'] = noise
#    fit_record['p_s']
#    fit_record['p_s_e']
    fit_record['v'] = np.array(fitted_vels,dtype=np.float32)
    fit_record['v_e'] = np.array(fitted_vels_e,dtype=np.float32)
    fit_record['w_l'] = np.array(fitted_width,dtype=np.float32)
    fit_record['w_l_e'] = np.array(fitted_width_e,dtype=np.float32) 
#    fit_record['w_s'] = 
#    fit_record['w_s_e'] =
#    fit_record['sd_l'] = 
#    fit_record['sd_s'] = 
#    fit_record['sd_phi'] = 
#    fit_record['x_qflg'] = 
#    fit_record['x_gflg'] = 
#    fit_record['x_p_l'] = 
#    fit_record['x_p_l_e'] = 
#    fit_record['x_p_s'] = 
#    fit_record['x_p_s_e'] = 
#    fit_record['x_v'] = 
#    fit_record['x_v_e'] = 
#    fit_record['x_w_l'] = 
#    fit_record['x_w_l_e'] = 
#    fit_record['x_w_s'] = 
#    fit_record['x_w_s_e'] = 
#    fit_record['phi0'] = 
#    fit_record['phi0_e'] = 
#    fit_record['elv'] = 
#    fit_record['elv_low'] = 
#    fit_record['elv_high'] = 
#    fit_record['x_sd_l'] = 
#    fit_record['x_sd_s'] = 
#    fit_record['x_sd_phi'] = 

    return fit_record


#Main program loop. Opens rawacf, fits it, writes fitted
#data to output filename.
def main(input_filename, output_filename):

    import numpy as np

    fitted_records = list()

    #first, read the file we want to fit
    records = read_file(input_filename)

    #next iterate through all records and fit them
    for record in records:
        fit_record = lmfit2(record)
        if isinstance(fit_record,dict):
            fitted_records.append(fit_record)
        break

    #finally, write the fitted records to a dmap file
    write_file(output_filename, fitted_records)


def lmfit2_mp(pos,record):
    return pos,lmfit2(record)


if __name__ == '__main__':
    import multiprocessing as mp

    import numpy as np
    pool = mp.Pool(10)

    input_filename = '20160303.1600.00.rkn.rawacf'
    output_filename = '20160303.1600.00.rkn.lmfit2'

    fitted_records = list()

    #first, read the file we want to fit
    records = read_file(input_filename)

    temp = [pool.apply_async(lmfit2_mp,args=(x,records[x])) for x in range(len(records))]

    fitted_records = [p.get() for p in temp]
    fitted_records.sort()
    fitted_records = [r[1] for r in fitted_records]

    write_file(output_filename, fitted_records)
