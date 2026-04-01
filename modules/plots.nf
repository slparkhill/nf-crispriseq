process plot_count_distributions {

   tag "${id}"
   label 'big_mem'
   cpus 1

   publishDir( 
      "${params.outputs}/plot-counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( umi_table )

   output:
   tuple val( id ), path( "*.{png,csv}" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd

   COUNTS_COLS = ["umi_count", "read_count"]
   
   df = (
      pd.read_csv(
         "${umi_table}", 
         sep='\\t',
         low_memory=False,
      )
   )
   cols = [col for col in COUNTS_COLS if col in df]
   if len(cols) == 1:
      cols *= 2

   scattergrid_args = {
      "df": df,
      "grid_columns": cols,
      "log": cols,
   }

   fig, axes = scattergrid(**scattergrid_args)
   figsaver(format="png")(
      fig=fig,
      name='count-hist',
      df=df,
   )

   """
   
}


process plot_counts_per_gene {

   tag "${id}"
   label 'big_mem'
   cpus 1

   publishDir( 
      "${params.outputs}/plot-counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( umi_table )

   output:
   tuple val( id ), path( "*.{png,csv}" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd
   
   df = (
      pd.read_csv(
         "${umi_table}", 
         sep='\\t',
      ) 
      .groupby(["genome_accession", "chr", "gene_id"])
   )

   bc_counts = (
      df[["umi_count", "read_count"]]
      .sum()
      .reset_index()
      .rename(columns={
         "umi_count": "umi_count_per_gene",
         "read_count": "read_count_per_gene"
      })
      .merge(
         df[["clone_barcode"]]
         .nunique()
         .reset_index()
      )
      .rename(columns={
         "clone_barcode": "clone_count_per_gene",
      })
   )

   cols = [
      "umi_count_per_gene", 
      "read_count_per_gene", 
      "clone_count_per_gene",
   ]

   fig, axes = scattergrid(
      bc_counts,
      grid_columns=cols,
      log=cols,
      grouping=["genome_accession"],
      aspect_ratio=1.25,
   )
   figsaver(format="png")(
      fig=fig,
      name='counts-per-gene',
      df=bc_counts,
   )

   """
   
}


process plot_counts_per_cell {

   tag "${id}"
   label 'big_mem'
   cpus 1

   publishDir( 
      "${params.outputs}/plot-counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( umi_table )

   output:
   tuple val( id ), path( "*.{png,csv}" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd

   GENE_ID = [
      "genome_accession", 
      "chr",
      "gene_id",
   ]

   df = pd.read_csv(
      "${umi_table}", 
      sep='\\t',
   )

   count_cols = ["umi_count", "read_count"]
   group_cols = GENE_ID + count_cols
   renamer = {
      col: f"{col}_per_cell" for col in group_cols
   }
   per_cell = (
      df
      .assign(
         gene_id=lambda x: x[GENE_ID[0]].str.cat(x[GENE_ID[1:]], sep=":"),
      )
      .groupby("cell_barcode")
      [group_cols]
      .agg({
         "gene_id": "nunique",
      } | {
         col: "sum" for col in count_cols
      })
      .rename(columns=renamer)
   )
   fig, axes = scattergrid(
      per_cell,
      grid_columns=list(renamer.values()),
      log=list(renamer.values()),
      #grouping=["gene_biotype"] if using_biotype else None,
      #aspect_ratio=1.25 if using_biotype else 1.,
   )
   figsaver(format="png")(
      fig=fig,
      name="counts-per-cell",
      df=per_cell.reset_index(),
   ) 

   """
   
}
