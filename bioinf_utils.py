import argparse
import os
from helpers.fastq import (
    parse_bounds,
    parse_fastq,
    is_gc_valid,
    is_length_valid,
    is_quality_valid
)
from helpers.bio_files import (
    read_fasta_lines,
    write_fasta_record,
    read_blast_file
)
def filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0.0):
    """
    Функция фильтрации. Сохраняет результат в папку filtered/
    """
    output_dir = "filtered"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_fastq)

    with open(output_path, 'w', encoding='utf-8') as out_f:
        for identifier, sequence, separator, quality in parse_fastq(input_fastq):
            if (is_gc_valid(sequence, gc_bounds) and
                    is_length_valid(sequence, length_bounds) and
                    is_quality_valid(quality, quality_threshold)):
                out_f.write(f"{identifier}\n{sequence}\n{separator}\n{quality}\n")

def convert_multiline_fasta_to_oneline(input_file_path: str, output_file_path: str = None):
    """
    Конвертирует многострочный FASTA файл в однострочный.
    Сохраняет результат в папку result/.
    """
    output_dir = "result"
    os.makedirs(output_dir, exist_ok=True)

    if output_file_path is None:
        base_name = os.path.basename(input_file_path)
        output_file_path = "converted_" + base_name

    output_file_path = os.path.join(output_dir, output_file_path)

    with open(output_file_path, 'w', encoding='utf-8') as outfile:
        current_header = None
        current_seq_parts = []

        for line in read_fasta_lines(input_file_path):
            if not line:
                continue

            if line.startswith('>'):
                if current_header is not None:
                    write_fasta_record(outfile, current_header, ''.join(current_seq_parts))
                current_header = line
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

        if current_header is not None:
            write_fasta_record(outfile, current_header, ''.join(current_seq_parts))


def parse_blast_output(blast_file_path: str, output_file_path: str = None):
    """
    Парсит результаты BLAST (текстовый формат). Для каждого запроса QUERY
     выбирает по первой строке из столбца `Description`. Сохраняет отсортированный по алфавиту результат в папку result/.
    """
    proteins = read_blast_file(blast_file_path)

    sorted_proteins = sorted(proteins)

    output_dir = "result"
    os.makedirs(output_dir, exist_ok=True)

    if output_file_path is None:
        base_name = os.path.basename(blast_file_path)
        output_file_path = "parsed_" + base_name

    output_file_path = os.path.join(output_dir, output_file_path)

    with open(output_file_path, 'w', encoding='utf-8') as out_f:
        for protein in sorted_proteins:
            out_f.write(protein + "\n")

    return sorted_proteins

def main():
    """
    Парсинг аргументов для всех утилит
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True, help="Доступные утилиты:")


    #FASTQ filtration
    parser_fastq = subparsers.add_parser("filter_fastq", help="Фильтрация FASTQ файла")
    parser_fastq.add_argument("-fq", required=True, help="Путь до входного FASTQ файла")
    parser_fastq.add_argument("-out", required=True, help="Имя выходного отфильтрованного файла. Отфильтрованные риды сохраняются \
        в папку filtered по имени `-out`")
    parser_fastq.add_argument("-gc", default=(0, 100), type=parse_bounds, help="Интервал GC-состава. По-умолчанию равен `(0, 100)`. \
        Если в аргумент передать одно число, то считается, что это верхняя граница. Примеры: `-gc (20, 80)` \
        - сохраняем только риды с GC составом от 20 до 80%%, `-gc 44.4` - сохраняем риды с GC составом меньше, чем 44.4%%")
    parser_fastq.add_argument("-len", default=(0, 2**32), type=parse_bounds, help="Интервал длины для фильтрации, всё аналогично '-gc', \
        но по-умолчанию равен `(0, 2**32)`")
    parser_fastq.add_argument("-qual", default=0.0, type=float, help="Пороговое значение среднего качества рида для фильтрации, \
        по-умолчанию равно `0` (шкала phred33). Риды со средним качеством по всем нуклеотидам ниже порогового отбрасываются.")

    #convert multiline FASTA to one line
    parser_fasta = subparsers.add_parser("convert_fasta", help="Конвертация многострочного FASTA файла в однострочный")
    parser_fasta.add_argument("-i", required=True, help="Путь до входного FASTA файла")
    parser_fasta.add_argument("-o", default=None, help="Имя выходного конвертированного файла (опционально)")

    #Parse BLAST output
    parser_blast = subparsers.add_parser("parse_blast", help="Парсинг и сортировка вывода BLAST")
    parser_blast.add_argument("-i", required=True, help="Путь до файла с результатами BLAST")
    parser_blast.add_argument("-o", default=None, help="Имя выходного отсортированного файла (опционально)")


    args = parser.parse_args()

    if args.command == "filter_fastq":
        # Call filtration function
        filter_fastq(
            input_fastq=args.fq,
            output_fastq=args.out,
            gc_bounds=args.gc,
            length_bounds=args.len,
            quality_threshold=args.qual
        )
    elif args.command == "convert_fasta":
        convert_multiline_fasta_to_oneline(args.i, args.o)

    elif args.command == "parse_blast":
        parse_blast_output(args.i, args.o)


if __name__ == "__main__":
    main()