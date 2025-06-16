; Reference: https://www.jianshu.com/p/eae362946ea9

; E[S_k]=\frac {k(k-1)}4，V_{ar}[S_k]=\frac {k(k-1)(2k+5)}{72} \quad 1≤k≤n
; S_k=\sum_{j=1}^kr_j \quad ,r_j=\begin{cases} 1&amp;X_j&gt;X_i\\ 0&amp;X_j≤X_i\\\end{cases} \quad (i=1,2,...,j;k=1,2,....n)
; UF_k=\frac {(S_k-E[S_k])}{\sqrt {V_{ar}[S_k]}}

function ukstd(sk, k)
; UF/UB  cacul procedure
begin
    uk=(sk-k*(k-1)/4)/sqrt(k*(k-1)*(2*k+5)/72.)
    return(uk)
end

function mktest(smp)
; 1-D Sample series only ; Length > 3
begin
    length=dimsizes(smp)
    skp=smp
    skp=0
    skn=skp
    uf=skp
    ub=skn
    smpr=smp(::-1)  ; Reversed sample series
    ; smpr(0)=smp(length-1)
    do numi = 2, length-1
        skp(numi)=skp(numi-1)+sum(where(smp(numi).gt.smp(:numi-1), 1, 0)) ;S_k of sample
        uf(numi)=ukstd(skp(numi), numi)

        skn(numi)=skn(numi-1)+sum(where(smpr(numi).gt.smpr(:numi-1), 1, 0)) ;S_k of reversed sample
        ub(numi)=ukstd(skn(numi), numi)
    end do

    ub=(/-ub(::-1)/) ;Reversed Back

    uf(:1)=-999.999
    uf(length-2:)=-999.999
    uf@_FillValue=-999.999
    ub(:1)=-999.999
    ub(length-2:)=-999.999
    ub@_FillValue=-999.999

    ufb=True
    ufb@uf=uf
    ufb@ub=ub

    return(ufb)
end

;; 最佳实践sample
begin
    sample=fspan(0, 5, 30)
    sample(15:)=fspan(5, 0, 15)
    ufb=mktest(sample)
    uf=ufb@uf
    ub=ufb@ub
    print(uf)
    print(ub)
end