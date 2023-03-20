using LaTeXTabulars
using LaTeXStrings
## drop rows with missing data
function na_omit!(data)
    #replace "missing" with "NA"
    df = replace(data,NaN    =>"NA")
    df = replace(df  ,missing=>"NA")
    df = replace(df  ,nothing=>"NA")
    df = replace(df  ,undef  =>"NA")

    sum(df .== "NA")
    #identify rows with "NA"
    na_rows     = sum(df.=="NA",dims=2) .== 0
    #enumerate rows
    total_rows  = collect(1:size(df,1))
    #identify rows that we want
    keep_rows   = convert.(Int64,na_rows) .* total_rows
    #reducing the number of rows
    filtred_rows=filter(x->x!=0, keep_rows)
    #replacement
    df2=df[filtred_rows,:]
    #@assert sum(df2 .== "NA") == 0 #check for "NA"
    return df2
end;