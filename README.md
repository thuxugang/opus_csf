# A series of scoring functions for ranking protein structural models

Our scoring functions are based on the native distributions (derived from the entire PDB) of coordinate components of specific atoms on selected residues of peptide segments of 5, 7, 9, and 11 residues in length, rather than using traditional Boltzmann formula. Therefore, the time complexity of each method is O(n), making them very suitable for applications at early stage of structural building. For more details, please see the reference papers.

## OPUS-CSF

OPUS-CSF is a fast and accurate scoring function that can be used to distinguish the native protein structures from their decoys based on the conformations of their ***main chains*** exclusively.

### Reference 
```bibtex
@article{xu2018opus,
  title={OPUS-CSF: A C-atom-based scoring function for ranking protein structural models},
  author={Xu, Gang and Ma, Tianqi and Zang, Tianwu and Wang, Qinghua and Ma, Jianpeng},
  journal={Protein Science},
  volume={27},
  number={1},
  pages={286--292},
  year={2018},
  publisher={Wiley Online Library}
}
```

## OPUS-DASF

OPUS-DASF is a fast and accurate scoring function that can be used to distinguish the native protein structures from their decoys based on the conformations of their ***side chains*** exclusively.


### Reference 
```bibtex
@article{xu2019opus,
  title={OPUS-Rota2: An Improved Fast and Accurate Side-chain Modeling Method},
  author={Xu, Gang and Ma, Tianqi and Du, Junqing and Wang, Qinghua and Ma, Jianpeng},
  journal={Journal of chemical theory and computation},
  year={2019},
  publisher={ACS Publications}
}
```

## OPUS-SSF

OPUS-SSF is a fast and accurate scoring function that can be used to distinguish the native protein structures from their decoys based on the conformations of their ***entire structures***.


### Reference 
```bibtex
@article{xu2019opus,
  title={OPUS-SSF: A side-chain-inclusive scoring function for ranking protein structural models},
  author={Xu, Gang and Ma, Tianqi and Wang, Qinghua and Ma, Jianpeng},
  journal={Protein Science},
  volume={28},
  number={6},
  pages={1157--1162},
  year={2019},
  publisher={Wiley Online Library}
}
```

## Performance

The results of OPUS-CSF, OPUS-SSF and OPUS-DASF on 5 decoy sets. The numbers of targets, with their native structures successfully recognized by each method, are listed in the table. The numbers in parentheses are the average Z-scores of the native structures. The larger the absolute value of Z-score, the better of our results.

||TOTAL|OPUS-CSF|OPUS-SSF|OPUS-DASF|
|:----:|:----:|:----:|:----:|:----:|
|3DRobot|200|189 (−4.86)	|186 (-5.24)|183 (-5.55)|
|Rosetta (3DR)|58|	51 (−3.83)	|53 (-3.98)	|52 (-3.95)|
|I-Tasser (3DR)|56|	36 (−3.47)	|38 (-3.81)	|38 (-3.40)|
|Rosetta|	58|	47 (−5.43)|	52 (-5.81)	|47 (-4.46)|
|I-Tasser|	56|	47 (−7.70)|	50 (-9.11)	|49 (-8.34)|


## Lookup Tables

To generate the four lookup tables (containing segments of 5, 7, 9, and 11 residues in length) in each method, we downloaded the entire PDB which contains 150,742 structures from `ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb` on May 10, 2019, and removed the structures in the test decoy sets.

The OPUS-CSF lookup tables are hosted on [Baidu Drive](https://pan.baidu.com/s/1OPDXv2y4C67KyN60Wk3X0w) with password `yqyy`. Also, they can be downloaded directly from [Here](http://ma-lab.rice.edu/dist/csf_lookup_tables.zip).

The OPUS-DASF lookup tables are hosted on [Baidu Drive](https://pan.baidu.com/s/1JmqS9T0kyRUqyeY-i3Pxdg ) with password `zw80`. Also, they can be downloaded directly from [Here](http://ma-lab.rice.edu/dist/dasf_lookup_tables.zip).

The OPUS-SSF lookup tables are hosted on [Baidu Drive](https://pan.baidu.com/s/1pqIJYjfq-kCi2f5gwL4SZQ) with password `9183`. Also, they can be downloaded directly from [Here](http://ma-lab.rice.edu/dist/ssf_lookup_tables.zip).

## Dependency

```
MongoDB v3.6.12
Python v3.6
pymongo v3.7.1
```

## Usage

1. Download the lookup tables you need and import them into your MongoDB. Then, create index. For example, you can use following scripts to import and index the csf lookup table with segments of 5 in length.

   ```
   mongoimport -d csf_db -c csf_5 --file your_path\csf_5.dat --type json
   use csf_db
   db.csf_5.ensureIndex({"key":1},{"unique":true})
   ```

2. Install pymongo.

   ```
   pip install pymongo
   ```

3. Set the parameters in *main.py* including `scoring_function` you would use and its corresponding `database_name` and `collections_name` in your MongoDB.

4. Run *main.py*.





