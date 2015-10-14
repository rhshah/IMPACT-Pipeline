Requirements
============

Please See template.conf file in the configuration folder.

:perl: `v5.20.2 <http://perl5.git.perl.org/perl.git/tag/2c93aff028f866699beb26e5e7504e531c31b284>`_
:python: `v2.7.8 <https://www.python.org/download/releases/2.7.8/>`_
:R: `v3.1.2 <http://cran.r-project.org/src/base/R-3/R-3.1.2.tar.gz>`_

Tools Used
----------
:Somatic SNV calling: `MuTect v1.1.4 <https://github.com/broadinstitute/mutect/tree/1.1.4>`_
:Somatic INDEL calling: `SomaticIndelDetector GATK v2.3-9 <http://www.broadinstitute.org/gatk/download>`_
:Somatic INDEL calling: `PINDEL v0.2.5a7 <https://github.com/genome/pindel/tree/v0.2.5a7>`_
:Somatic Structural Variant Framework: `IMPACT-SV v1.0.1 <https://github.com/rhshah/IMPACT-SV/tree/1.0.1>`_

Inside the config file
----------------------

There are three sections:

+-----------+-----------+-----------+
| Section 1 | Section 2 | Section 3 |
+===========+===========+===========+
| Locations | Parameters| Versions  |
+-----------+-----------+-----------+

Inside Location Here are the things that need to be set:

:ZCAT: Location of the ``zcat`` program on linux 
:TMPDIR: Set the temporary directory for all tools please set somthing other then ``/tmp``
:JAVA_1_6: Set JAVA version 1.6
:JAVA_1_7: Set JAVA version 1.7
:GATK_SomaticIndel: Path to GATK somatic indel detector (GATK version 2.3-9)
:GATK: Path to GATK (GATK version 3.3.0)
:Reference: Path to fasta referece file to be used (GRCh37)
:Refseq: Path to refgene file to be used
:PICARD: Path to picard tools (Picard version 1.19)

Set the parameters to different file/folders/values required by the IMPACT pipeline

Inside the version there are version that are being used for each tool. This is just for consistency in reports.

   
   