miRA a micro RNA identification tool
=========

Software and source code for the paper:
*Conservation-independent identification of novel miRNAs*
M. Evers, A. Dueck, G. Meister and J. C. Engelmann

### How to install:

If you are running a recent 64bit Linux, or OSX Lion or later you can grab the miRA binary here:

- [Linux-64bit](bin/miRA-Linux-64bit)
- [OSX-64bit](bin/miRA-OSX-64bit)

You will still need gnuplot and latex installed for your system. 
You will also need the Varna binary: [Varna](VARNAv3-91.jar)
Please place it into the same Folder the miRA binary is in. 

#### Compiling from source

For an easy setup simply download the 
latest bundeled release archive: [miRA-1.1.1.tar.gz](miRA-1.1.1.tar.gz)


unpack it, using for example
```sh
tar -xvf miRA-1.1.1.tar.gz
```
Make sure your system supplies the following dependecies for miRA:

- a c compiler supporting the c99 standard
- a java virtual machine version 1.6+ (optional)
- a recent version of gnuplot (optional)
- a recent version of latex (optional)

**NOTE:** miRA will work without the optional dependencies but will skip some reporting features (creating plots etc.) if they are not available.

Compile it for your system with:
```sh
cd miRA-1.0.5
./configure
make
```

Optionally run the unit tests on your system with: 
```sh
make test
```
to check for correct behavior.




### How to use:

The simplest and most common way to run miRA is to run the full Suite using the command:
```sh
./miRA full -c <configuration file> <input SAM file> <input FASTA file> <output directory>
```

You can test miRA with sample data provided in [./example/](example):
```sh
./miRA full -c example/sample_configuration.config example/sample_reads.sam example/sample_sequence.fasta example/sample_output/
```

You can also run only parts of miRA, it is seperated in 3 parts with distinct calls for each one:

| Algorithm        | Description           | Command  |
| ------------- |:-------------:| :--------------|
| Clustering     | generates a list of main expression contigs based on alignment data| cluster |
| Folding    | fold rna sequences and calculate secondary structure information      |   fold |
| Coverage Testing | coverage based verification and reporting of micro rna candidates     |    coverage | 

For additional help and usage information run:
```sh
./miRA <command> -h
```
where <command\> is either "cluster" "fold" or "coverage"

###Results:

After running miRA all result files will be created in the specified output directory. Depending on the configuration and the available external programs the following files will be created:

- a full pdf report for every microRNA candidate (requires latex)
- final_candidates.bed, a file containing location and properties of all candidates in the bed file format.
- final_candidaes.json, a file containing location and properties of all candidates in the json file format.




