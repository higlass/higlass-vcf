# HiGlass VCF Track

> Track for viewing VCF files.

[![HiGlass](https://img.shields.io/badge/higlass-👍-red.svg?colorB=0f5d92)](http://higlass.io)
[![Build Status](https://img.shields.io/travis/higlass/higlass-pileup-track/master.svg?colorB=0f5d92)](https://travis-ci.org/higlass/higlass-pileup-track)

<img src="/teaser.png?raw=true" width="600" />

**Note**: This is the source code for the pileup only! You might want to check out the following repositories as well:

- HiGlass viewer: https://github.com/higlass/higlass
- HiGlass server: https://github.com/higlass/higlass-server
- HiGlass docker: https://github.com/higlass/higlass-docker

## Installation

```
npm install higlass-vcf
```

## Usage

The live scripts can be found at:

- https://unpkg.com/higlass-vcf@v0.1.2/dist/higlass-vcf.min.js

### Client

1. Make sure you load this track prior to `hglib.js`. For example:

```
<script src="https://unpkg.com/higlass-vcf@v0.1.2/dist/higlass-vcf.min.js"></script>
<script src="hglib.js"></script>
<script>
  ...
</script>
```

2. Now, configure the track in your view config and be happy!

```
{
  "editable": true,
  "trackSourceServers": [
    "http://higlass.io/api/v1"
  ],
  "exportViewUrl": "/api/v1/viewconfs",
  "views": [
    {
      "initialXDomain": [
        0,
        100000
      ],
      "tracks": {
        "top": [
          {
            "type": "vcf",
            "height": 180,
            "uid": "FylkvVBTSumoJ959HT4-5A",
            "data": {
              "type": "vcf",
              "url": "https://pkerp.s3.amazonaws.com/public/HG002_SVs_Tier1_v0.6.vcf.gz",
              "chromSizesUrl": "https://resgen.io/api/v1/chrom-sizes/?id=cNE4StljSAK9lK3amECl-A"
            },
          }
        ]
      }
    }
  ]
}
```

## Support

For questions, please either open an issue or ask on the HiGlass Slack channel at http://bit.ly/higlass-slack

## Development

### Installation

```bash
$ git clone https://github.com/higlass/higlass-vcf-track && higlass-vcf-track
$ npm install
```

### Commands

**Developmental server**: `npm start`
**Production build**: `npm run build`
