
<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/abdrakhimov1/Biosaur">
    <img src="logo.png" alt="Logo" width="506" height="222">
  </a>

  <h3 align="center">Biosaur</h3>

  <p align="center">
    Modern software for data analysis
    <br />
    <a href="https://github.com/abdrakhimov1/Biosaur"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/abdrakhimov1/Biosaur">View Demo</a>
    ·
    <a href="https://github.com/abdrakhimov1/Biosaur/issues">Report Bug</a>
    ·
    <a href="https://github.com/abdrakhimov1/Biosaur/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Installation](#installation)
* [Usage](#usage)
* [Targeted Mode](#targeted-mode)
* [Roadmap](#roadmap)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)



<!-- ABOUT THE PROJECT -->
## About The Project


`Biosaur`: open source peptide MS feature detector.

`Biosaur` provides the opportunity to work with:
* Data captured in negative mode
* Data containing information about ion mobility
* Also `biosaur`reports the correlation map

`Biosaur` algorithm allows users to get all the functionality of standard isotope detecting tool with the additional ability to analyze ion mobility data from devices of different types (such as FIAMS TimsTOF)

### Built With
Biosaurus was developed using
* [python3](https://www.python.org/)
* [numpy](https://numpy.org/)
* [pandas](https://pandas.pydata.org/)



<!-- GETTING STARTED -->
## Getting Started

`Biosaur` is a console utility that is easy to install and configure on your personal computer or computing cluster.

### Installation

There are several options to install `Biosaur`.

* Easy way: you can use 
```sh
pip3 install biosaur 
```
which inastall stable version of `biosaur` on your computer. 

* If you want to get latest actual version of `biosaur` you shold use next algorithm:

1. Clone the repo
```sh
git clone https://github.com/abdrakhimov1/Biosaur.git
```
3. Enter the `biosaur` directory
```sh
cd Biosaur
```
4. Install `biosaur`:
```sh
pip3 insatall .
```

<!-- USAGE EXAMPLES -->
## Usage

Biosaur is quite easy to use. To start your first bioasur search use command:
```sh
biosaur YOUR_FILE.mzML
```
This command will start standart biosaur search with default parameters.
If you need to specify parameters use `biosaur --help` to identify the required parameter.

Special attention to TIMS TOF data.

First of all, the .d files should be converted to mzML format using msconvert with option '--combineIonMobilitySpectra'.

**Please, do not use option `--filter "scanSumming"`! The latter is often required for MS/MS data analysis but breaks MS1 feature detection.**

The best way to deal with it is to use `--combineIonMobilitySpectra` with `--filter "msLevel 1"` to create an individual mzML file for Biosaur-only analysis. At the current moment, TIMS TOF data has enormous size of files, as well as a huge amount of peaks, so it is highly recommended to use Biosaur `--min_intensity` option to reduce complexity of the analysis. For example, using `--min_intensity 1000` option requires ~10 Gb of RAM memory and 20 mins of processing time on average PC (Intel i7-3930K CPU) when applied to a complex sample dataset containing 8000 MS1 spectra (`200ng_HeLa_50cm_120min_100ms from PXD010012` on the ProteomeXchange). The same data with `--min_intensity 800` filter requires 40 minutes of processing. The analyis of similar data for Orbitrap HF with no ion mobility info and no restrictions on `--min_intensity` takes ~5-10 min. In general, increasing `--min_intensity` reduces Biosaur analysis time and RAM consumption in non-linear way, but at the same time decreases sensitivity of feature detection.

## Targeted Mode

Biosaur has targeted mode, in which it matches the results of identification of MS/MS spectra to the peptide features. To activate it, the MS/MS search results in pepXML or mzID format are required. Biosaur will take into account MS/MS search results during feature detection workflow.
If you want to activate biosaur targeted mode, add a keyword `--pxfp`  and provide path to the results of the MS/MS search engine.

Current Biosaur version supports X!Tandem, IdentiPy, MSFragger, Comet search outputs in pepXML formats, as well as MSGF+ output in mzID format.

Example:
```sh
biosaur YOUR_FILE.mzML --pxfp YOUR_SEARCH_ENGINE_RESULT.pep.xml
```
The output of Biosaur will contain a column with the MS/MS scans IDs and the corresponding peptide features, as in the standart mode of the Biosaur.

<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/abdrakhimov1/Biosaur/issues) for a list of proposed features (and known issues).



<!-- CONTRIBUTING -->
## Contributing

We are open and **welcome** various collaborations with representatives of the international community so we are ready to discuss any improvements to the biosaur. Any contributions you make are **greatly appreciated**. To help us with `biosaur` improvment follow next steps.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request
6. Contact us for discussion



<!-- LICENSE -->
## License

Distributed under the Apache 2.0 License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

* Abdrakhimov Daniil - [@abdrakhimov1](https://t.me/abdrakhimov1) - dan.abdrakhimov@yandex.ru
* Ivanov Mark - [@markmipt](https://t.me/markmipt) - markmipt@gmail.com

Project Link: [https://github.com/abdrakhimov1/Biosaur](https://github.com/abdrakhimov1/Biosaur)


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/abdrakhimov1/Biosaur
[contributors-url]: https://github.com/abdrakhimov1/Biosaur/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/abdrakhimov1/Biosaur
[forks-url]: https://github.com/abdrakhimov1/Biosaur/network/members
[stars-shield]: https://img.shields.io/github/stars/abdrakhimov1/Biosaur
[stars-url]: https://github.com/abdrakhimov1/Biosaur/stargazers
[issues-shield]: https://img.shields.io/github/issues/abdrakhimov1/Biosaur
[issues-url]: https://github.com/abdrakhimov1/Biosaur/issues
[license-shield]: https://img.shields.io/github/license/abdrakhimov1/Biosaur
[license-url]: https://www.apache.org/licenses/LICENSE-2.0


