process plot_UMI_distributions {

   tag "${id}"
   label 'big_mem'

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
   
   df = pd.read_csv(
      "${umi_table}", 
      sep='\\t',
      low_memory=False,
   )

   fig, axes = scattergrid(
      df,
      grid_columns=["umi_count", "bulk_read_count"],
      log=["umi_count", "bulk_read_count"],
      grouping=["gene_biotype"],
      aspect_ratio=1.25,
   )

   figsaver(format="png")(
      fig=fig,
      name='umi-hist',
      df=df,
   )

   """
   
}
