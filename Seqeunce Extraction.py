#!/usr/bin/env python3
"""Fetch MEROPS substrate annotations for protease UniProt IDs and export CSV/HTML results.

This single script combines the earlier workflow into one reusable tool:
- accepts either a plain list of UniProt IDs or a UniProt TSV export
- fetches MEROPS cross-references from UniProt
- fetches MEROPS substrate pages and parses substrate rows
- fetches full protein sequences from UniProt when available
- writes CSV and HTML outputs
"""

import argparse
import csv
import html
import json
import re
import ssl
import time
import urllib.request
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, List, Optional, Tuple

UNIPROT_URL = 'https://rest.uniprot.org/uniprotkb/{}.json'
MEROPS_SUBSTRATE_URL = 'https://www.ebi.ac.uk/merops/cgi-bin/substrates?id={}'

OUTPUT_HEADERS = [
    'Input_UniProt',
    'Protease_Sequence',
    'MEROPS_ID',
    'Substrate',
    'Uniprot',
    'Residue range',
    'Cleavage Site',
    'Cleavage type',
    'Evidence',
    'P4',
    'P3',
    'P2',
    'P1',
    "P1'",
    "P2'",
    "P3'",
    "P4'",
    'Substrate_Sequence',
    'Reference',
    'CutDB',
    'MERNUM',
]


class TableCell:
    def __init__(self) -> None:
        self.text = ''
        self.links: List[Dict[str, str]] = []


class MeropsTableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.rows: List[List[TableCell]] = []
        self.current_row: Optional[List[TableCell]] = None
        self.current_cell: Optional[TableCell] = None
        self.in_cell = False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'tr':
            self.current_row = []
        elif tag in ('td', 'th') and self.current_row is not None:
            self.current_cell = TableCell()
            self.current_row.append(self.current_cell)
            self.in_cell = True
        elif tag == 'a' and self.current_cell is not None:
            href = attrs.get('href', '')
            if href:
                self.current_cell.links.append({'href': href})

    def handle_data(self, data):
        if self.in_cell and self.current_cell is not None:
            self.current_cell.text += data

    def handle_endtag(self, tag):
        if tag in ('td', 'th') and self.in_cell:
            self.in_cell = False
            if self.current_cell is not None:
                self.current_cell.text = self.current_cell.text.strip()
                self.current_cell = None
        elif tag == 'tr' and self.current_row is not None:
            if any(cell.text or cell.links for cell in self.current_row):
                self.rows.append(self.current_row)
            self.current_row = None


def fetch_url(url: str) -> str:
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(url, headers={'User-Agent': 'python-urllib/3'})
    with urllib.request.urlopen(req, context=ctx, timeout=30) as response:
        return response.read().decode('utf-8', errors='ignore')


def read_input_ids(path: str, column: Optional[str] = None) -> List[str]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f'Input file not found: {path}')

    with p.open('r', encoding='utf-8', errors='ignore', newline='') as f:
        first = f.readline()
        if not first:
            return []

        if '\t' in first or ',' in first:
            f.seek(0)
            try:
                reader = csv.DictReader(f, delimiter='\t' if '\t' in first else ',')
            except Exception:
                reader = csv.DictReader(f)

            fieldnames = [name.strip() for name in (reader.fieldnames or [])]
            if column and column in fieldnames:
                return [row.get(column, '').strip() for row in reader if row.get(column, '').strip()]
            if 'Entry' in fieldnames:
                return [row.get('Entry', '').strip() for row in reader if row.get('Entry', '').strip()]
            if 'UniProt' in fieldnames:
                return [row.get('UniProt', '').strip() for row in reader if row.get('UniProt', '').strip()]
            if 'UniProt_ID' in fieldnames:
                return [row.get('UniProt_ID', '').strip() for row in reader if row.get('UniProt_ID', '').strip()]
            raise ValueError(f'Could not find an ID column in {path}. Available columns: {fieldnames}')

        if first.strip():
            return [first.strip()] + [line.strip() for line in f if line.strip()]
        return []


def extract_uniprot_tsv_fields(path: str) -> List[Dict[str, str]]:
    p = Path(path)
    with p.open('r', encoding='utf-8', errors='ignore', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = []
        for row in reader:
            entry = (row.get('Entry') or row.get('Entry ID') or '').strip()
            if entry:
                rows.append({
                    'Entry': entry,
                    'Sequence': (row.get('Sequence') or '').strip(),
                    'Entry_Name': (row.get('Entry Name') or '').strip(),
                    'Protein_names': (row.get('Protein names') or '').strip(),
                    'Organism': (row.get('Organism') or '').strip(),
                })
        return rows


def fetch_merops_ids_from_uniprot(uniprot_id: str) -> List[str]:
    url = UNIPROT_URL.format(uniprot_id)
    try:
        data = json.loads(fetch_url(url))
        return [ref['id'] for ref in data.get('dbReferences', []) if ref.get('type') == 'MEROPS']
    except Exception:
        return []


def fetch_protein_sequence(uniprot_id: str) -> str:
    url = UNIPROT_URL.format(uniprot_id)
    try:
        data = json.loads(fetch_url(url))
        sequence_obj = data.get('sequence', {})
        if isinstance(sequence_obj, dict):
            return sequence_obj.get('value', '')
        return ''
    except Exception:
        return ''


def parse_merops_substrates(html_text: str) -> List[Dict[str, str]]:
    parser = MeropsTableParser()
    parser.feed(html_text)
    if not parser.rows:
        return []

    headers = [cell.text for cell in parser.rows[0]]
    uniprot_index = next((i for i, header in enumerate(headers) if header.lower() == 'uniprot'), None)
    results: List[Dict[str, str]] = []
    for row in parser.rows[1:]:
        values = [cell.text for cell in row]
        row_dict = {headers[i]: values[i] if i < len(values) else '' for i in range(len(headers))}
        if uniprot_index is not None and uniprot_index < len(row):
            uniprot_cell = row[uniprot_index]
            if not row_dict.get('Uniprot') and uniprot_cell.links:
                for link in uniprot_cell.links:
                    if 'SpAcc=' in link['href']:
                        match = re.search(r'SpAcc=([A-Za-z0-9_-]+)', link['href'])
                        if match:
                            row_dict['Uniprot'] = match.group(1)
                            break
        elif len(row) > 1:
            uniprot_cell = row[1]
            if not row_dict.get('Uniprot') and uniprot_cell.links:
                for link in uniprot_cell.links:
                    if 'SpAcc=' in link['href']:
                        match = re.search(r'SpAcc=([A-Za-z0-9_-]+)', link['href'])
                        if match:
                            row_dict['Uniprot'] = match.group(1)
                            break
        results.append(row_dict)
    return results


def fetch_substrates_for_merops_id(merops_id: str) -> List[Dict[str, str]]:
    html_text = fetch_url(MEROPS_SUBSTRATE_URL.format(merops_id))
    return parse_merops_substrates(html_text)


def write_csv(path: str, rows: List[Dict[str, str]]) -> None:
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_HEADERS)
        writer.writeheader()
        if rows:
            writer.writerows(rows)


def write_html(path: str, rows: List[Dict[str, str]]) -> None:
    with open(path, 'w', encoding='utf-8') as f:
        f.write('<!doctype html>\n')
        f.write('<html lang="en">\n')
        f.write('<head>\n')
        f.write('  <meta charset="utf-8">\n')
        f.write('  <title>Protease Substrate Results</title>\n')
        f.write('  <style>body{font-family:Arial,sans-serif;margin:20px;}table{border-collapse:collapse;width:100%;}th,td{border:1px solid #ccc;padding:6px 8px;vertical-align:top;}th{background:#f2f2f2;text-align:left;}tr:nth-child(even) td{background:#fbfbfb;}</style>\n')
        f.write('</head>\n')
        f.write('<body>\n')
        f.write('<h1>Protease Substrate Results</h1>\n')
        f.write('  <table>\n')
        f.write('    <thead><tr>\n')
        for header in OUTPUT_HEADERS:
            f.write(f'      <th>{html.escape(header)}</th>\n')
        f.write('    </tr></thead>\n')
        f.write('    <tbody>\n')
        for row in rows:
            f.write('      <tr>\n')
            for header in OUTPUT_HEADERS:
                value = row.get(header, '')
                f.write(f'        <td>{html.escape(str(value), quote=True)}</td>\n')
            f.write('      </tr>\n')
        f.write('    </tbody>\n')
        f.write('  </table>\n')
        f.write('</body>\n')
        f.write('</html>\n')


def build_rows(uniprot_ids: List[str], delay: float) -> List[Dict[str, str]]:
    sequence_cache: Dict[str, str] = {}
    output_rows: List[Dict[str, str]] = []

    for uniprot_id in uniprot_ids:
        print('Processing', uniprot_id)
        protease_sequence = sequence_cache.get(uniprot_id, '')
        if not protease_sequence:
            protease_sequence = fetch_protein_sequence(uniprot_id)
            sequence_cache[uniprot_id] = protease_sequence

        merops_ids = fetch_merops_ids_from_uniprot(uniprot_id)
        if not merops_ids:
            print('  No MEROPS cross-reference found for', uniprot_id)
            continue

        for merops_id in merops_ids:
            time.sleep(delay)
            try:
                substrates = fetch_substrates_for_merops_id(merops_id)
            except Exception as exc:
                print(f'  Failed to fetch substrates for {merops_id}: {exc}')
                continue
            if not substrates:
                print('  No substrates found for', merops_id)
                continue

            for substrate in substrates:
                substrate_uniprot = substrate.get('Uniprot', '')
                substrate_sequence = ''
                if substrate_uniprot:
                    substrate_sequence = sequence_cache.get(substrate_uniprot, '')
                    if not substrate_sequence:
                        time.sleep(delay)
                        substrate_sequence = fetch_protein_sequence(substrate_uniprot)
                        sequence_cache[substrate_uniprot] = substrate_sequence

                output_rows.append({
                    'Input_UniProt': uniprot_id,
                    'Protease_Sequence': protease_sequence,
                    'MEROPS_ID': merops_id,
                    'Substrate': substrate.get('Substrate', ''),
                    'Uniprot': substrate.get('Uniprot', ''),
                    'Residue range': substrate.get('Residue range', ''),
                    'Cleavage Site': substrate.get('Cleavage Site', ''),
                    'Cleavage type': substrate.get('Cleavage type', ''),
                    'Evidence': substrate.get('Evidence', ''),
                    'P4': substrate.get('P4', ''),
                    'P3': substrate.get('P3', ''),
                    'P2': substrate.get('P2', ''),
                    'P1': substrate.get('P1', ''),
                    "P1'": substrate.get("P1'", ''),
                    "P2'": substrate.get("P2'", ''),
                    "P3'": substrate.get("P3'", ''),
                    "P4'": substrate.get("P4'", ''),
                    'Substrate_Sequence': substrate_sequence,
                    'Reference': substrate.get('Reference', ''),
                    'CutDB': substrate.get('CutDB', ''),
                    'MERNUM': substrate.get('MERNUM', ''),
                })
    return output_rows


def main() -> int:
    parser = argparse.ArgumentParser(description='Acquire protease substrate annotations from UniProt and MEROPS.')
    parser.add_argument('input', help='Input file: plain IDs, CSV/TSV with an ID column, or a UniProt TSV export.')
    parser.add_argument('output', nargs='?', default='protease_substrates.csv', help='Output CSV file path.')
    parser.add_argument('--html-output', help='Optional HTML output path.')
    parser.add_argument('--column', '-c', default=None, help='Column name to read UniProt IDs from when the input is tabular.')
    parser.add_argument('--delay', type=float, default=1.0, help='Seconds to wait between web requests.')
    args = parser.parse_args()

    try:
        input_rows = extract_uniprot_tsv_fields(args.input)
        if input_rows:
            uniprot_ids = [row['Entry'] for row in input_rows]
            print(f'Loaded {len(uniprot_ids)} UniProt entries from TSV export: {args.input}')
        else:
            uniprot_ids = read_input_ids(args.input, args.column)
            print(f'Loaded {len(uniprot_ids)} IDs from: {args.input}')
    except Exception as exc:
        print('Error reading input:', exc)
        return 2

    if not uniprot_ids:
        print('No UniProt IDs found in input.')
        return 1

    output_rows = build_rows(uniprot_ids, args.delay)

    write_csv(args.output, output_rows)
    html_output = args.html_output or (args.output[:-4] + '.html' if args.output.lower().endswith('.csv') else args.output + '.html')
    write_html(html_output, output_rows)

    print('Wrote', len(output_rows), 'rows to', args.output)
    print('Wrote', len(output_rows), 'rows to', html_output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
