import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy.optimize import curve_fit
R=8.314 # gas constant [kg m^2 / K mol s^2]

def fit_isotherm(st_data,con_data,min_pt, max_pt,n,σ_sol,T,isotherm):
    l1=min_pt # lower limit index for low concentration (pre CMC) data
    l2=max_pt# upper limit index for low concnetration (pre CMC) dat
    # split points into isotherm fitting region and after CMC region
    st1=st_data[l1:l2] #N/m      
    con1=con_data[l1:l2]
    st2=st_data[l2:len(st_data)]
    con2=con_data[l2:len(con_data)]
    def langmuir(surf_con,a,b): # Equate Eq 10 and 12 from Estoe and Dalton (2000) Advnaces in Colloid and Interface Science, 85, 103.
        return σ_sol+ a*np.log(1-((b*surf_con)/(1+b*surf_con)))
    def rsquared(x, y):
        """ Return R^2 where x and y are array-like."""
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
        return r_value**2
    if isotherm == 'gibs':
        pp=np.polyfit(np.log(con1), st1, 1,cov=True) # fit low concnetration data surface tension vs ln(concentration) with straight line
        em1=np.sqrt(pp[1][0][0])
        eb1=np.sqrt(pp[1][1][1])
        x1=np.linspace(round(np.log(min(con_data))),round(np.log(max(con_data))),100)
        sigmadat=pp[0][0]*x1+pp[0][1]
        Gamma=-1*(pp[0][0])/(n*R*T) #gibbs isotherm to find surface excess using derivative of linear fit above
        Gamma_error=1*(eb1)/(n*R*T) #gibbs isotherm to find surface excess using derivative of linear fit above
        pl=np.polyfit(con2, st2, 1,cov=True) #fit high concentration data with straight line
        em2=np.sqrt(pl[1][0][0])
        eb2=np.sqrt(pl[1][1][1])
        y=pl[0][0]*x1+pl[0][1] # fit line
        intersect=(pp[0][1]-pl[0][1])/(pl[0][0]-pp[0][0]) #find intersection of low and high concentration fits == CMC 
        cmc_error=np.exp(intersect)*np.sqrt(((eb2**2+eb1**2)/(pl[0][1]-pp[0][1])**2)+((em1**2+em2**2)/(pl[0][0]-pp[0][0])**2)) # error on intersetion of two stratight lines
        r2=rsquared(np.log(con1),st1)
        full_con=np.exp(x1)
        plt.plot(x1,sigmadat)
        plt.scatter(np.log(con_data),st_data)
        plt.plot(x1,y)
        plt.ylabel('Surface Tension [N/m]')
        plt.xlabel('ln(concentration / mol m^-3)')
        plt.legend(['Pre CMC linear fit','Beyond CMC linear fit'])
        plt.ylim([min(st_data[l1:len(st_data)])-0.1*min(st_data[l1:len(st_data)]),max(st_data[l1:len(st_data)])*1.1])
        plt.show()
        print('CMC-fits=',np.exp(intersect),'+/-',cmc_error,'mM','Gammamax=',Gamma,'+/-', Gamma_error,'mol/m^2')
        dat={'CMC (mM)': np.exp(intersect),
            'CMCerror (mM)':cmc_error,
            'G_max (mol/m^2)':Gamma,
            'G_max error (mol/m^2)':[Gamma_error],
        }
        fitdat={'concentration (mM)':full_con,
            'isotherm fit (N/m)':sigmadat,
            'plateau fit (N/m)':y}
        df=pd.DataFrame(dat)
        df.to_csv('Gibs_Isotherm_Parameters_Output.csv')
        df2=pd.DataFrame(fitdat)
        df2.to_csv('Gibs_Isotherm_Fit_Output.csv')
    elif isotherm == 'langmuir':
        popt,pcov=curve_fit(langmuir,con1,st1)
        variance = np.diagonal(pcov)
        standard_error=np.sqrt(variance)
        G_max=popt[0]/(n*R*T)
        G_maxx_error=standard_error[0]/(n*R*T)
        a=popt[1]
        a_error=standard_error[1]
        #x_fit=np.linspace(0,max(con_data),100)
        x_fit=np.linspace(np.log10(min(con_data)),np.log10(max(con_data)),500) # updated to have better spacing of points over orders of magnitude
        y_fit=langmuir(x_fit,popt[0],popt[1])
        pl=np.polyfit(con2, st2, 1,cov=True) #fit high concentration data with straight line
        em2=np.sqrt(pl[1][0][0])
        eb2=np.sqrt(pl[1][1][1])
        y=pl[0][0]*x_fit+pl[0][1] # fit line
        dif=np.array(y_fit-y) # determine CMC using minimum difference between two region fitting lines
        idx,diff=min(enumerate(dif), key=lambda x: abs(x[1]-0))
        cmc=x_fit[idx]
        cmc_error=(x_fit[2]-x_fit[1])/2 # error on CMC set to 1/2 the spacing of the points for for the fit lines
        plt.plot(x_fit,y_fit)
        plt.plot(x_fit,y)
        plt.scatter(con_data,st_data)
        plt.ylabel('Surface Tension [N/m]')
        plt.xlabel('Concentration  (mol m^-3)')
        plt.legend(['Pre CMC fit','Beyond CMC fit','Measurements'])
        plt.ylim([min(st_data[l1:len(st_data)])-0.1*min(st_data[l1:len(st_data)]),max(st_data[l1:len(st_data)])*1.1])
        plt.show()
        print('CMC=',cmc,'+/-',cmc_error,'mM G_max=',G_max,'+/-',G_maxx_error,'mol/m^2 a=',popt[1],'+/-',a_error,'m^3/mol')
        dat={'CMC (mM)': cmc,
            'CMCerror (mM)':cmc_error,
            'G_max (mol/m^2)':G_max,
            'G_max error (mol/m^2)':G_maxx_error,
            'a (m^3/mol)':a,
            'a error (m^3/mol)':[a_error],
        }
        fitdat={'concentration (mM)':x_fit,
            'isotherm fit (N/m)':y_fit,
            'plateaus fit (N/m)':y}
        df=pd.DataFrame(dat)
        df.to_csv('Langmuir_Isotherm_Parameter_Output.csv')
        df2=pd.DataFrame(fitdat)
        df2.to_csv('Langmuir_Isotherm_Fit_Output.csv')
