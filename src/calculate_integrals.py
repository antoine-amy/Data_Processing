def calculate_integrals_simple(data, inf, sup):
    nn=len(data)
    AREAs= []
    wsize=data.shape[1]
    print(wsize)
    hsize=int(wsize/2)
    for i in range(nn):
        dtl, dtr, tfit, tlim, tll = -15, 110, 350, 600, 30
        bl=np.mean(data[i][hsize+2*dtl:hsize+dtl])
        wf=data[i]-bl
        maxx=np.max(wf)
        area=np.sum(wf)
        
        ww, hh=8, 12
        dled=wf[ww:]-wf[:-ww]
        listpeaks,_=find_peaks(dled, height=hh,distance=20)
        peakpos=listpeaks[(listpeaks<sup) & (listpeaks>inf)]
        AREAs.append(area)
    data = pd.DataFrame(columns=['area'])
    data['area'] = AREAs
    return data