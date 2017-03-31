# NGI-RNAseq with AWS

## Introduction

Amazon Web Services (AWS) are cloud-based compute systems, suitable for large-scale bioinformatics analysis. They can be an attractive way to run scalable analyses without needing to pay large up-front expenses for hardware and maintenance.

NGI-RNAseq is built with NextFlow, which has a numer of AWS integrations built in. However, AWS is vast and almost infinitely configurable, and as such there are several different ways in which you can use it for your analyses. Which approach you use depends on your needs - the frequency with which you run the pipeline, the scale of data you need to analysis and your familiarity with running in the cloud.

In this document, we describe the steps taken when running NGI-RNAseq on AWS hardware using three different methods. These have been written as a walk-through based on our experience at the time of writing. Note that we have used AWS based in Ireland (`eu-west-1`) and that the availability of different hardware can vary according to location.

If you find any problems with this documentation, please let us know or better still, submit a pull-request with a fix.

## Table of Contents

1. [AWS Concepts](#aws-concepts)
    * Basic introduction to some AWS terms that you'll need to know.
2. [Fully Manual Single Node Setup](#1-fully-manual-single-node-setup)
    * The most basic way to run the pipeline. Involves creating a server on AWS and running NextFlow without any clever parallelisation stuff.
    * Good when you're just getting started, or if you have a small amount of data.
3. NextFlow Integration - Single Node _(coming soon)_
4. NextFlow Integration - Elastic Cluster _(coming soon)_

## AWS Concepts
The people behind AWS love acronyms. There is a huge amount of jargon to get your head around when starting with this. We've tried to keep this documentation as simple as possible, but it helps to understand a few basics:

### Instance Types
Broadly speaking, everything in this document will use compute resources that are part of the AWS EC2 - that is, the Amazon Web Services Elastic Cloud 2. These are on-demand servers that you hire on a use-by-time basis. Each one you run is called an _instance_. You can choose different instance types according to your requirements as each has different hardware specifications. See [Amazon EC2 Instance Types](https://aws.amazon.com/ec2/instance-types/) for a full list.

### Storage types
As well as an instance to run computation work on, we'll also need some file storage. This is organised separately from the instance and is then mounted to it. There are three types which are worth knowing about:

1. EFS _(Elastic File System)_
    * This is the most expensive filesystem. It can scales to whatever size of data you need, without any need to predefine a requirement. You just pay for what you use. A blank EFS filesystem is free to keep alive.
    * EFS filesystems can be mounted on multiple different instances.
2. EBS _(Elastic Block Store)_
    * This filesystem is cheaper than EFS. You can choose what type of storage to use (eg. solid-state or hard disk) for varying prices.
    * Unlike EFS, you need to know how much storage will be required in advance. If you go over this then everything will break.
    * EBS volumes can only be mounted on a single instance.
3. S3 _(Simple Storage Service)_
    * Cheap long term storage.
    * It's good for archiving data to once you're finished, but no good for running any analysis on.

In this documentation, we use EFS for running all of our analysis. There are a few benefits it gives - if the instance dies, we don't lose our data. It's faster and we don't need to worry about how much space we're going to need. The downside is that it's expensive, so you should clear it out as soon as you're done.

### Regions
AWS resources are based in data centres, which are based in regions. Typically, all data and hardware is restricted to a single region. These examples are written based on experience with Ireland. They should work in US-West with some tweaking. Other regions may have varying hardware available, so approach with care.

## 1. Fully Manual Single Node Setup
Simplest to setup, just using a single instance.

1. Log into your AWS console
    * Select the Ireland region (`eu-west-1`)
    * This can all be done either with the web interface or the command line client. Here we will describe the web interface.
2. Launch an EFS File system
    1. Go to the _Storage: EFS_ Dashboard
    2. Click _Create Storage_
    3. Leave everything as default and click through the pages to create the filesystem.
        * Note the name of the security group that is created
    4. Once created, expand the EFS details and copy the _DNS name_ to your clipboard for later.
3. Set up a new security group
    1. Go to the EC2 Dashboard
    2. Select _Security Groups_ under the _Network & Security_ section in the left navigation
    3. Click _Create Security Group_
    4. Set the name and description to `SSH` (or whatever you like)
    5. Click _Add Rule_
    6. Select `SSH` as _Type_ and `Anywhere` as _Source_ (or a more limited IP range if you prefer)
    7. Click _Create_
4. Launch an Amazon Machine Instance (AMI)
    1. Go to the EC2 Dashboard
    2. Click _Launch Instance_
    3. Select the _Community AMIs_ tab and search for the NGI-RNAseq AMI that we have created
        * AMI id: `ami-9a5c79e9`
        * This image contains all of the required software for the NGI-RNAseq pipeline and is maintained by us for easy use by the community.
    4. Click Select
    5. Select the `m4.4xlarge` Instance Type
        * If running with STAR and Human data, we need around 34GB RAM. This is the cheapest instance type with that memory available.
        * If aligning with HISAT2 or using a smaller genome then you can use smaller instances, such as `m4.2xlarge` or `t2.2xlarge` (`t2` is less efficient but cheaper).
    6. Click _Configure Instance_
    7. Optional: Click _Request Spot Instances_
        * Enter a maximum bid value. Choosing something close to or above the on-demand price will mean that you rarely if ever get kicked off, but still get reduced costs most of the time.
        * Using spot pricing typically means that you save around 87% of the list price. See the [AWS Spot Bid Advisor](https://aws.amazon.com/ec2/spot/bid-advisor/) for more information.
        * For `m4.4xlarge`, $0.95 should be sufficient at time of writing.
    8. Select one of the subnets and note the name.
    9. Expand _Advanced Details_ at the bottom.
        * Under _User data_ add the following as _as text_:
        ```yaml
        # cloud-config
        package_upgrade: true
        packages:
        - nfs-utils
        runcmd:
        - mkdir -p /mnt/efs
        - chown ec2-user:ec2-user /mnt/efs
        - echo "SUBNET-ID.YOUR-DNS-NAME-HERE:/ /mnt/efs nfs4 nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2 0 0" >> /etc/fstab
        - mount -a -t nfs4
        - sudo chown ec2-user /mnt/efs
        ```
        * Replace `SUBNET-ID` and `YOUR-DNS-NAME-HERE` with the subnet ID from above and the DNS name that you copied earlier.
        * For example, you may copy in something like: `eu-west-1c.fs-hwjer3hs.efs.eu-west-1.amazonaws.com`
    10. Click through the tabs until you get to _Configure Security Group_
    11. Use _Select an existing security group_ and choose the two that you created earlier
        * One to access the EFS, one to SSH into the instance.
    12. Click _Review and Launch_, then _Launch_
    13. Create a new Key Pair if you don't already have one
        * This is used to access the instance by `ssh`
    14. Click _Request Spot Instances_ :tada:
5. Access the new instance
    1. Once the spot instance has been given, go to _Instances_ on the EC2 dashboard.
    2. Copy the _IPv4 Public IP_ from the new instance details tab
    3. Load up a console on your local computer (eg. _Terminal_)
    4. Once the _Instance State_ and _Status Checks_ are complete, then SSH into the new instance
        * `ssh -i /path/to/key.pem ec2-user@YOUR-IP`
        * Replace `YOUR-IP` with the _IPv4 Public IP_ that you copied above.
    5. If you manage to log in successfully, pat yourself on the back and have a cookie. :cookie:
6. Get Nextflow running with the pipeline
    1. Go to the `efs` directory and work there: `cd /mnt/efs`
    2. Download the required references
    3. Download the input data
        * You can try this with our test data:
        ```
        curl https://export.uppmax.uu.se/b2013064/test-data/ngi-rna_test_set.tar.bz2 > ngi-rna_test_set.tar.bz2
        tar xvjf ngi-rna_test_set.tar.bz2
        ```
    4. Run the pipeline!
        * If you've downloaded the test data, this command should work:
        ```
        nextflow run SciLifeLab/NGI-RNAseq -profile base --gtf ngi-rna_test_set/genes.gtf --bed12 ngi-rna_test_set/genes.bed --star_index ngi-rna_test_set/star/ --reads "ngi-rna_test_set/*.fastq.gz"
        ```
6. Retrieve your finished results to your local computer
    1. Use SCP on your local machine to pull your data back from the EFS file system:
    ```
    scp -r -i /path/to/key.pem ec2-user@YOUR-IP:/mnt/efs/results .
    ```
7. Stop everything and clean up
    1. Wipe the EFS filesystem so that it doesn't cost any money
        * On the AWS instance, delete all files on EFS:
        ```
        cd ~
        rm -rf /mnt/efs/*
        ```
        * Now it's empty, it should be free to keep it around (so no need to delete it)
    2. Go to the EC2 Dashboard
    3. Select the running instance and click _Actions_ > _Instance State_ > _Terminate_
        * **If you leave this running, things will get expensive really quickly**
    4. Cancel the spot request under _Spot Requests_

