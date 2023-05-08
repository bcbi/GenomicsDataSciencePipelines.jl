using ConfParser
using CSV
using DataFrames
using Dates
using StatsBase
using UnicodePlots

# Unique values in df.col, sorted by number of occurrences
function sorted_counts(df, col)
	sort(combine(groupby(df, col), nrow), :nrow, rev=true)
end

"""
# Find data dictionary entry for a column
# data_dictionary("Conditions") # Example
"""
function data_dictionary(x)
	ddf[ddf[!,"Variable Name"] .== uppercase(x),:]
end

# Read in options
conf = ConfParse("./config.ini")
parse_conf!(conf)
data_dir = retrieve(conf, "local", "data_dir")
prior_file = retrieve(conf, "local", "prior_file")
variant_file = retrieve(conf, "local", "variant_file")

work_dir = retrieve(conf, "local", "work_dir")
dict_file = retrieve(conf, "local", "dict_file")
accession_file = retrieve(conf, "local", "accession_file")

# Read in data
vdf = CSV.File(joinpath(data_dir, variant_file),
	       dateformat="mm/dd/yyyy",
	       types=Dict(:idl_zipcode=>String, :test_site_zip=>String) ) |> DataFrame
ddf = CSV.read(joinpath(work_dir, dict_file), DataFrame, select=1:6)
gdf = CSV.read(joinpath(work_dir, accession_file), DataFrame)

# Delete columns where all values are the same
for df in [vdf, gdf]
	bad_columns = names(df)[ df |> eachcol .|> allequal |> findall ]
	select!(df, Not( bad_columns ))
end

# Set proper datatype for true/false survey data
for col in names(vdf)
	all_values = unique(vdf[!,col])
	if all_values == [0,1]
		vdf[!,col] = Bool.(vdf[!,col])
	elseif all_values == ["No", "Yes"]
		d = Dict(["No", "Yes"] .=> [false,true])
		vdf[!,col] = vdf[!,col] .|> x->d[x]
	end
end

# For consistency between data sets, use all caps for column names
rename!(vdf, names(vdf) .=> uppercase.(names(vdf)))
rename!(gdf, names(gdf) .=> uppercase.(names(gdf)))

# Add ACCESSION_NUMBER to GISAID data matching the COVID-19 accession IDs
gdf.ACCESSION_NUMBER = let
	m = match.(r"^hCoV-19/USA/RI-RISHL-(.*)/20([0-9]{2})$",gdf.STRAIN)
	x = convert(Vector{Union{Nothing,String}}, fill(nothing,nrow(gdf)))
	for i in 1:nrow(gdf)
		a=m[i]
		isnothing(a) && continue
		x[i] = "RI-" * a[2] * "-" * a[1]
	end
	x;
end

# Join datasets on matching ACCESSION_NUMBER
jdf = innerjoin(vdf, gdf, on=:ACCESSION_NUMBER, matchmissing=:notequal)

# One row of gdf seems to have AGE and SEX interchanged.
# However, it doesn't appear in the innerjoin, so it is ignored.
