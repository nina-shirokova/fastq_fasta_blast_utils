import argparse
import os
from helpers.fastq import (
    parse_bounds,
    parse_fastq,
    is_gc_valid,
    is_length_valid,
    is_quality_valid
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


def main():
    """
    Парсинг аргументов
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-fq", required=True, help="Путь до входного fastq файла")
    parser.add_argument("-out", required=True, help="Имя выходного отфильтрованного файла. Отфильтрованные риды сохраняются \
        в папку filtered по имени `-out`")
    parser.add_argument("-gc", default=(0, 100), type=parse_bounds, help="Интервал GC-состава. По-умолчанию равен `(0, 100)`. \
        Если в аргумент передать одно число, то считается, что это верхняя граница. Примеры: `-gc (20, 80)` \
        - сохраняем только риды с GC составом от 20 до 80%%, `-gc 44.4` - сохраняем риды с GC составом меньше, чем 44.4%%")
    parser.add_argument("-len", default=(0, 2**32), type=parse_bounds, help="Интервал длины для фильтрации, всё аналогично '-gc', \
        но по-умолчанию равен `(0, 2**32)`")
    parser.add_argument("-qual", default=0.0, type=float, help="Пороговое значение среднего качества рида для фильтрации, \
        по-умолчанию равно `0` (шкала phred33). Риды со средним качеством по всем нуклеотидам ниже порогового отбрасываются.")

    args = parser.parse_args()

    #Call filtration function
    filter_fastq(
        input_fastq=args.fq,
        output_fastq=args.out,
        gc_bounds=args.gc,
        length_bounds=args.len,
        quality_threshold=args.qual
    )

if __name__ == "__main__":
    main()
print("OK")