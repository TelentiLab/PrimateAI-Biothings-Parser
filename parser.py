import logging
import os

FILE_NOT_FOUND_ERROR = 'Cannot find input file: {}'  # error message constant

# configure logger
logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s', level=logging.INFO)
_logger = logging.getLogger('biothings_parser')

# change following parameters accordingly
SOURCE_NAME = 'primate_ai'  # source name that appears in the api response
FILENAME = 'sample_data.tsv'  # name of the file to read
DELIMITER = '\t'  # the delimiter that separates each field


def _inspect_file(filename: str) -> int:
    i = 0
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def load_data(data_folder: str):
    """
    Load data from a specified file path. Parse each line into a dictionary according to the schema.
    Then process each dict by normalizing data format, remove null fields (optional).
    Append each dict into final result using its id.

    :param data_folder: the path(folder) where the data file is stored
    :return: a generator that yields data.
    """
    input_file = os.path.join(data_folder, FILENAME)
    # raise an error if file not found
    if not os.path.exists(input_file):
        _logger.error(FILE_NOT_FOUND_ERROR.format(input_file))
        raise FileExistsError(FILE_NOT_FOUND_ERROR.format(input_file))

    file_lines = _inspect_file(input_file)  # get total lines so that we can indicate progress in next step

    with open(input_file, 'r') as file:
        _logger.info(f'start reading file: {FILENAME}')
        count = 0
        skipped = []
        for line in file:
            count += 1
            _logger.info(f'reading line {count} ({(count / file_lines * 100):.2f}%)')  # format to use 2 decimals

            if line.startswith('#') or line.strip() == '':
                skipped.append(line)
                continue  # skip commented/empty lines

            try:
                (chrom, pos, ref, alt, ref_aa, alt_aa, strand_1pos_0neg, trinucleotide_context, ucsc_gene,
                 exac_coverage, primate_dl_score) = line.strip().split(DELIMITER)  # unpack according to schema
            except ValueError:
                _logger.error(f'failed to unpack line {count}: {line}')
                _logger.error(f'got: {line.strip().split(DELIMITER)}')
                skipped.append(line)
                continue  # skip error line

            try:  # parse each field if necessary (format, enforce datatype etc.)
                chrom = chrom.replace('chr', '')
                pos = int(pos)
                strand_1pos_0neg = bool(strand_1pos_0neg)
                exac_coverage = float(exac_coverage)
                primate_dl_score = float(primate_dl_score)
            except ValueError as e:
                _logger.error(f'failed to cast type for line {count}: {e}')
                skipped.append(line)
                continue  # skip error line

            _id = f'chr{chrom}:g.{pos}_{pos}'  # define id

            variant = {
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'ref_aa': ref_aa,
                'alt_aa': alt_aa,
                'pos_strand': strand_1pos_0neg,
                'trinucleotide_context': trinucleotide_context,
                'ucsc_gene': ucsc_gene,
                'exac_coverage': exac_coverage,
                'primate_dl_score': primate_dl_score,
            }

            yield {  # commit an entry by yielding
                "_id": _id,
                SOURCE_NAME: variant
            }
        _logger.info(f'parse completed, {len(skipped)}/{file_lines} lines skipped.')
        for x in skipped:
            _logger.info(f'skipped line: {x.strip()}')
