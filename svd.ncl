;-------------------------------------------;
;                                           ;
;             SVD Analysis                  ;
;                                           ;
;   This file repacked a new function to do ;
; 3d svd analysis using svdcov more easily. ;
;                                           ;
;   Input: 2 matrix & how many mods need    ;
;                                           ;
;   Output: HOMOGENEOUS matrixes & ts       ;
;           varience                        ;
;                                           ;
;-------------------------------------------;



function svdfunc(lfield,rfield,nsvd)

    begin
        timestmp0=systemfunc("date +%s")
        nyr=dimsizes(lfield(:,0,0))

        lmod=ind(.not.ismissing(ndtooned(lfield(0,:,:))))
        lmod_nn=dimsizes(ndtooned(lmod))
        ; printVarSummary(lmod)

        rmod=ind(.not.ismissing(ndtooned(rfield(0,:,:))))
        rmod_nn=dimsizes(ndtooned(rmod))
        ; printVarSummary(rmod)

        lfield_1d=ndtooned(lfield)
        ind_lf1=ind(.not.ismissing(lfield_1d))
        ; printVarSummary(ind_lf1)
        x_1d=lfield_1d(ind_lf1)
        ; printVarSummary(x_1d)
        x_2d=onedtond(x_1d,(/nyr,lmod_nn/))
        x_2d!0="year"
        x_2d!1="grid"
        x_2d:=x_2d(grid|:,year|:)

        rfield_1d=ndtooned(rfield)
        ind_rf1=ind(.not.ismissing(rfield_1d))
        ; printVarSummary(ind_rf1)
        y_1d=rfield_1d(ind_rf1)
        y_2d=onedtond(y_1d,(/nyr,rmod_nn/))
        y_2d!0="year"
        y_2d!1="grid"
        y_2d:=y_2d(grid|:,year|:)

        homlft=new((/nsvd,lmod_nn/),float)
        hetlft=new((/nsvd,lmod_nn/),float)
        homrgt=new((/nsvd,rmod_nn/),float)
        hetrgt=new((/nsvd,rmod_nn/),float)
        svdrst=svdcov(x_2d, y_2d, nsvd, homlft, hetlft, homrgt, hetrgt)
        pcvar=svdrst;ceil(svdrst*10)/10.0
        svd_ts=new((/nsvd*2,nyr/), float)
        svd_ts!0="mods"
        svd_ts&mods=ispan(0, (nsvd*2-1), 1)
        svd_ts(:(nsvd-1),:)=onedtond(svdrst@ak, (/nsvd,nyr/))
        svd_ts(nsvd:,:)=onedtond(svdrst@bk, (/nsvd,nyr/))

        homleft_n1d=ndtooned(lfield(:nsvd-1,:,:))
        lnsvdmod=ind(.not.ismissing(homleft_n1d))
        homleft_n1d=homleft_n1d@_FillValue
        homleft_n1d(lnsvdmod)=ndtooned(homlft)
        lphase=onedtond(homleft_n1d, dimsizes(lfield(:nsvd-1,:,:)))

        homrgt_n1d=ndtooned(rfield(:nsvd-1,:,:))
        rnsvdmod=ind(.not.ismissing(homrgt_n1d))
        homrgt_n1d=homrgt_n1d@_FillValue
        homrgt_n1d(rnsvdmod)=ndtooned(homrgt)
        rphase=onedtond(homrgt_n1d, dimsizes(rfield(:nsvd-1,:,:)))

        timestmp1=systemfunc("date +%s")
        print("Function runned in "+(timestmp1-timestmp0)+" seconds.")

        return([/lphase,rphase,svd_ts,pcvar/])
    end