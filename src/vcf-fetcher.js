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

     // console.log('chromInfo', chromInfo);

    const tbiIndexed = new TabixIndexedFile({
      filehandle: new RemoteFile(this.dataConfig.url),
      tbiFilehandle: new RemoteFile(`${this.dataConfig.url}.tbi`) // can also be csiFilehandle
    });

    const headerText = await tbiIndexed.getHeader();
    // console.log('headerText:', headerText);

    const TILE_SIZE = 1024;
    const tbiVCFParser = new VCF({ header: headerText })

    // console.log('tbiVCFParser', tbiVCFParser);

    const retVal = {
      tile_size: TILE_SIZE,
      bins_per_dimension: TILE_SIZE,
      max_zoom: Math.ceil(
        Math.log(chromInfo.totalLength / TILE_SIZE) / Math.log(2)
      ),
      max_width: chromInfo.totalLength,
      min_pos: [0],
      max_pos: [chromInfo.totalLength]
    };

    this.tbiIndexed =  tbiIndexed
    this.chromInfo = chromInfo;
    this.tbiVCFParser = tbiVCFParser;
    this.tsInfo =  retVal

    console.log('retVal:', retVal);
    callback(retVal);
  }


  async tile(z, x) {
    const MAX_TILE_WIDTH = 200000;
    const tsInfo= this.tsInfo;
    console.log('tsInfo:', tsInfo);
    const tileWidth = +tsInfo.max_width / 2 ** +z;
    const fetched =  [];

    if (tileWidth > MAX_TILE_WIDTH) {
      // this.errorTextText('Zoomed out too far for this track. Zoomin further to see reads');
      return [];
    }

    // get the bounds of the tile
    let minX = tsInfo.min_pos[0] + x * tileWidth;
    const maxX = tsInfo.min_pos[0] + (x + 1) * tileWidth;

    const chromInfo = this.chromInfo;
    const results = [];

    const { chromLengths, cumPositions } = chromInfo;

    for (let i = 0; i < cumPositions.length; i++) {
      const chromName = cumPositions[i].chr;
      const chromStart = cumPositions[i].pos;

      const chromEnd = cumPositions[i].pos + chromLengths[chromName];

      if (chromStart <= minX && minX < chromEnd) {
        // start of the visible region is within this chromosome

        if (maxX > chromEnd) {
          const  startPos = minX - chromStart
          const  endPos = chromEnd - chromStart

          await this.tbiIndexed.getLines(chromName,
            startPos, endPos, line =>
            fetched.push(this.tbiVCFParser.parseLine(line)),
          );
          console.log('sp  1', startPos, endPos);
          // the visible region extends beyond the end of this chromosome
          // fetch from the start until the end of the chromosome
          // continue onto the next chromosome
          minX = chromEnd;
        } else {
          const endPos = Math.ceil(maxX - chromStart);
          const startPos = Math.floor(minX - chromStart);
          // the end of the region is within this chromosome

          await this.tbiIndexed.getLines(chromName,
            startPos, endPos, line =>
            fetched.push(this.tbiVCFParser.parseLine(line)),
          );
          console.log('sp 2', startPos, endPos);
          // end the loop because we've retrieved the last chromosome
          break;
        }
      }
    }

    // flatten the array of promises so that it looks like we're
    // getting one long list of value
    console.log('fetched:', fetched);
    return fetched;
  }

  fetchTilesDebounced(receivedTiles, tileIds) {
    const tiles = {};

    const validTileIds = [];
    const tilePromises = [];

    for (const tileId of tileIds) {
      const parts = tileId.split('.');
      const z = parseInt(parts[0], 10);
      const x = parseInt(parts[1], 10);

      if (Number.isNaN(x) || Number.isNaN(z)) {
        console.warn('Invalid tile zoom or position:', z, x);
        continue;
      }


      // validTileIds.push(tileId);
      // tilePromises.push(tile(uid, z, x));

      this.tile(z, x);
    }

    // return Promise.all(tilePromises).then(values => {
    //   for (let i = 0; i < values.length; i++) {
    //     const validTileId = validTileIds[i];
    //     tiles[validTileId] = values[i];
    //     tiles[validTileId].tilePositionId = validTileId;
    //   }

    //   return tiles;
    // });
  }
}

export default VCFDataFetcher;
