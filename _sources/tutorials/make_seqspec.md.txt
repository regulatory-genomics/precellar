# Write a seqspec file for your custom datasets

## Introduction to seqspec
Precellar relies on a seqspec file to understand the structure of your datasets. For details about seqspec, please refer to the [seqspec documentation](https://github.com/pachterlab/seqspec).

We have provided some seqspec templates for common datasets in the [`seqspec/templates` directory](https://github.com/pachterlab/precellar/tree/main/seqspec/templates). If you have a dataset that is not yet supported, you can check whether it is provided [here](https://www.sina.bio/seqspec-builder/assays.html).


## Notes on writing seqspec files

### Format
The seqspec file is a yaml file that adheres to the [seqspec schema](https://github.com/pachterlab/seqspec/blob/main/seqspec/schema/seqspec.schema.json). 

### Customizing seqspec files
Sometimes the library structure varies from the default templates. For example, barcodes length and location can be different in different Cel-seq2 datasets. In this case, you can customize the seqspec file to match your dataset.
