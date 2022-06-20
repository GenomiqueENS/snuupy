import pandas as pd
import argparse
from tqdm import tqdm

#parse args
parser = argparse.ArgumentParser()
parser.add_argument('-i', metavar='mismatch results', dest='MISMATCH_RESULT', nargs='+', required=True,
                        help='getMismatch output')
parser.add_argument( "-o", dest='OUTPUT_FEATHER',  help="nanopore read id with barcode and umi; feather format")
parser.add_argument( "--ED-barcode", dest='MAX_BARCODE_ED', type= int, help="MAX BARCODE ED")
parser.add_argument( "--ED-UMI", dest='MAX_UMI_ED', type= int, help="MAX_UMI_ED")



def barcodeAssignment(mismatchResult, MAX_BARCODE_ED, MAX_UMI_ED):
    """
    assign barcode for each Nanopore read; based on mismatch results
    """

    mismatchResult.drop_duplicates(['name','qseqid'], inplace=True)

    mismatchResult['barcodeUmiMismatch'] = mismatchResult['barcodeUmiMismatch'].astype(int)
    mismatchResult['barcodeMismatch'] = mismatchResult['barcodeMismatch'].astype(int)
    mismatchResult['umiMismatch'] = mismatchResult['umiMismatch'].astype(int)

    mismatchResult['hitContent'] = mismatchResult['qseqid'] + ',' + mismatchResult['barcodeUmiMismatch'].astype(str) + ','\
        + mismatchResult['barcodeMismatch'].astype(str) + ',' + mismatchResult['umiMismatch'].astype(str) + ';'
    
    mismatchResult['allHitContents'] = mismatchResult.groupby('name')['hitContent'].transform('sum')
    mismatchResult['hitCounts'] = mismatchResult.groupby('name')['hitContent'].transform('count')

    mismatchResult = mismatchResult.loc[(mismatchResult['barcodeMismatch'] <= MAX_BARCODE_ED) & \
            (mismatchResult['umiMismatch'] <= MAX_UMI_ED)]
    mismatchResult.sort_values(['name','barcodeUmiMismatch','barcodeMismatch','umiMismatch'], inplace=True)
    mismatchResult.drop_duplicates('name', inplace=True)
    
    mismatchResult.reset_index(drop=True, inplace=True)
    return mismatchResult

def save_feather(mismatchResult, OUTPUT_FEATHER):
    mismatchResult.to_feather(OUTPUT_FEATHER)


def main ():
    args = parser.parse_args()
    df_Assignment = pd.DataFrame()
    flag = True
    print('merging feathre files ...') 
    for fl in tqdm(args.MISMATCH_RESULT, total = len(args.MISMATCH_RESULT)):
        if flag:
            df_tmp = pd.read_feather(fl)
            df_Assignment = barcodeAssignment(df_tmp, args.MAX_BARCODE_ED, args.MAX_UMI_ED)
            flag = False
            continue
        df_tmp = pd.read_feather(fl)
        df_tmp = barcodeAssignment(df_tmp, args.MAX_BARCODE_ED, args.MAX_UMI_ED)
        df_Assignment = df_Assignment.append(df_tmp, ignore_index=True)
    print('merging files success!')
    print('barcode assignment ...')
    save_feather(df_Assignment, args.OUTPUT_FEATHER)
    df_Assignment = barcodeAssignment(df_Assignment, args.MAX_BARCODE_ED, args.MAX_UMI_ED)
    save_feather(df_Assignment, str(args.OUTPUT_FEATHER[:-8]+"_reduced.feather"))
    print('done!')
if __name__ == "__main__":
    main()
