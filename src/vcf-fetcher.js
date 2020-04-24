import { TabixIndexedFile } from '@gmod/tabix';
import VCF from '@gmod/vcf';
import {RemoteFile} from 'generic-filehandle';

import { tsvParseRows } from 'd3-dsv';

import { bisector } from 'd3-array';
import { text } from 'd3-request';

/////////////////////////////////////////////////
/// ChromInfo
/////////////////////////////////////////////////

const chromInfoBisector = bisector(d => d.pos).left;
const segmentFromBisector = bisector(d => d.from).right;

const chrToAbs = (chrom, chromPos, chromInfo) =>
  chromInfo.chrPositions[chrom].pos + chromPos;

const absToChr = (absPosition, chromInfo) => {
  if (!chromInfo || !chromInfo.cumPositions || !chromInfo.cumPositions.length) {
    return null;
  }

  let insertPoint = chromInfoBisector(chromInfo.cumPositions, absPosition);
  const lastChr = chromInfo.cumPositions[chromInfo.cumPositions.length - 1].chr;
  const lastLength = chromInfo.chromLengths[lastChr];

  insertPoint -= insertPoint > 0 && 1;

  let chrPosition = Math.floor(
    absPosition - chromInfo.cumPositions[insertPoint].pos
  );
  let offset = 0;

  if (chrPosition < 0) {
    // before the start of the genome
    offset = chrPosition - 1;
    chrPosition = 1;
  }

  if (
    insertPoint === chromInfo.cumPositions.length - 1 &&
    chrPosition > lastLength
  ) {
    // beyond the last chromosome
    offset = chrPosition - lastLength;
    chrPosition = lastLength;
  }

  return [
    chromInfo.cumPositions[insertPoint].chr,
    chrPosition,
    offset,
    insertPoint
  ];
};

function natcmp(xRow, yRow) {
  const x = xRow[0];
  const y = yRow[0];

  if (x.indexOf('_') >= 0) {
    const xParts = x.split('_');
    if (y.indexOf('_') >= 0) {
      // chr_1 vs chr_2
      const yParts = y.split('_');

      return natcmp(xParts[1], yParts[1]);
    }

    // chr_1 vs chr1
    // chr1 comes first
    return 1;
  }

  if (y.indexOf('_') >= 0) {
    // chr1 vs chr_1
    // y comes second
    return -1;
  }

  const xParts = [];
  const yParts = [];

  for (const part of x.match(/(\d+|[^\d]+)/g)) {
    xParts.push(isNaN(part) ? part.toLowerCase() : +part);
  }

  for (const part of y.match(/(\d+|[^\d]+)/g)) {
    xParts.push(isNaN(part) ? part.toLowerCase() : +part);
  }

  // order of these parameters is purposefully reverse how they should be
  // ordered
  for (const key of ['m', 'y', 'x']) {
    if (y.toLowerCase().includes(key)) return -1;
    if (x.toLowerCase().includes(key)) return 1;
  }

  if (xParts < yParts) {
    return -1;
  } else if (yParts > xParts) {
    return 1;
  } else {
    return 0;
  }

  return 0;
}

function parseChromsizesRows(data) {
  const cumValues = [];
  const chromLengths = {};
  const chrPositions = {};

  let totalLength = 0;

  for (let i = 0; i < data.length; i++) {
    const length = Number(data[i][1]);
    totalLength += length;

    const newValue = {
      id: i,
      chr: data[i][0],
      pos: totalLength - length
    };

    cumValues.push(newValue);
    chrPositions[newValue.chr] = newValue;
    chromLengths[data[i][0]] = length;
  }

  return {
    cumPositions: cumValues,
    chrPositions,
    totalLength,
    chromLengths
  };
}

function ChromosomeInfo(filepath, success) {
  const ret = {};

  ret.absToChr = absPos => (ret.chrPositions ? absToChr(absPos, ret) : null);

  ret.chrToAbs = ([chrName, chrPos] = []) =>
    ret.chrPositions ? chrToAbs(chrName, chrPos, ret) : null;

  return text(filepath, (error, chrInfoText) => {
    if (error) {
      // console.warn('Chromosome info not found at:', filepath);
      if (success) success(null);
    } else {
      const data = tsvParseRows(chrInfoText);
      const chromInfo = parseChromsizesRows(data);

      Object.keys(chromInfo).forEach(key => {
        ret[key] = chromInfo[key];
      });
      if (success) success(ret);
    }
  });
}

/////////////////////////////////////////////////////
/// End Chrominfo
/////////////////////////////////////////////////////


class VCFDataFetcher {
  constructor(dataConfig, HGC) {
    this.dataConfig = dataConfig;
    this.uid = HGC.libraries.slugid.nice();
    this.isServerFetcher = !(dataConfig.type && dataConfig.type === 'bam');
  }

  async tilesetInfo(callback) {
    const [chromSizesUrl, vcfUrl ] = [
      this.dataConfig.chromSizesUrl,
      this.dataConfig.url
     ];

     const chromInfo = await new Promise((resolve, reject) => {
       ChromosomeInfo(chromSizesUrl, resolve)
     });

     console.log('chromInfo', chromInfo);

    const tbiIndexed = new TabixIndexedFile({
      filehandle: new RemoteFile(this.dataConfig.url),
      tbiFilehandle: new RemoteFile(`${this.dataConfig.url}.tbi`) // can also be csiFilehandle
    });

    const headerText = await tbiIndexed.getHeader();
    console.log('headerText:', headerText);

    // const promises = chromSizesUrl
    //   ? [bamHeaders[bamUrl], chromSizes[chromSizesUrl]]
    //   : [bamHeaders[bamUrl]];

    // return Promise.all(promises).then(values => {
    //   const TILE_SIZE = 1024;
    //   let chromInfo = null;

    //   if (values.length > 1) {
    //     // we've passed in a chromInfo file
    //     // eslint-disable-next-line prefer-destructuring
    //     chromInfo = values[1];
    //   } else {
    //     // no chromInfo provided so we have to take it
    //     // from the bam file index
    //     const chroms = [];
    //     for (const { refName, length } of bamFiles[bamUrl].indexToChr) {
    //       chroms.push([refName, length]);
    //     }

    //     chroms.sort(natcmp);

    //     chromInfo = parseChromsizesRows(chroms);
    //   }

    //   chromInfos[chromSizesUrl] = chromInfo;

    //   const retVal = {
    //     tile_size: TILE_SIZE,
    //     bins_per_dimension: TILE_SIZE,
    //     max_zoom: Math.ceil(
    //       Math.log(chromInfo.totalLength / TILE_SIZE) / Math.log(2)
    //     ),
    //     max_width: chromInfo.totalLength,
    //     min_pos: [0],
    //     max_pos: [chromInfo.totalLength]
    //   };

    //   return retVal;
    // });

    return retVal = {
      tile_size: TILE_SIZE,
      bins_per_dimension: TILE_SIZE,
      max_zoom: Math.ceil(
        Math.log(chromInfo.totalLength / TILE_SIZE) / Math.log(2)
      ),
      max_width: chromInfo.totalLength,
      min_pos: [0],
      max_pos: [chromInfo.totalLength]
    };
  }

  fetchTilesDebounced(receivedTiles, tileIds) {
    // this.track.updateLoadingText();

    // this.worker.then(tileFunctions => {
    //   if (this.isServerFetcher) {
    //     tileFunctions
    //       .serverFetchTilesDebounced(this.uid, tileIds)
    //       .then(receivedTiles);
    //   } else {
    //     tileFunctions
    //       .fetchTilesDebounced(this.uid, tileIds)
    //       .then(receivedTiles);
    //   }
    // });
  }
}

export default VCFDataFetcher;
