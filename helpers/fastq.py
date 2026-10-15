import argparse

def parse_bounds(value):
    """
    Вспомогательная функция для парсинга границ GC-состава или длины
    из аргументов командной строки.
    """
    if value.startswith('(') and value.endswith(')'):
        value = value[1:-1]

    parts = value.split(',')
    if len(parts) == 1:
        try:
            return float(parts[0]) if '.' in parts[0] else int(parts[0])
        except ValueError:
            raise argparse.ArgumentTypeError(f"Неверный формат числа: {value}")
    elif len(parts) == 2:
        try:
            low = float(parts[0].strip()) if '.' in parts[0] else int(parts[0].strip())
            high = float(parts[1].strip()) if '.' in parts[1] else int(parts[1].strip())
            return (low, high)
        except ValueError:
            raise argparse.ArgumentTypeError(f"Неверный формат интервала: {value}")
    raise argparse.ArgumentTypeError(f"Неверный формат: {value}")

def parse_fastq(file_path):
    """
    Читает fastq-файл построчно блоками по 4 строки.
    Возвращает (yield) кортеж из 4 переменных для каждого рида.
    """
    with open(file_path, 'r', encoding='utf-8') as file:
        while True:
            identifier = file.readline()
            if not identifier:
                break
            sequence = file.readline()
            separator = file.readline()
            quality = file.readline()

            yield identifier.strip(), sequence.strip(), separator.strip(), quality.strip()


def is_gc_valid(sequence, gc_bounds):
    """
    Проверяет, попадает ли GC-состав рида в заданный интервал.
    """
    if not sequence:
        return False
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    gc_percent = (gc_count / len(sequence)) * 100

    if isinstance(gc_bounds, (int, float)):
        return gc_percent <= gc_bounds
    else:
        return gc_bounds[0] <= gc_percent <= gc_bounds[1]

def is_length_valid(sequence, length_bounds):
    """
    Проверяет, попадает ли длина рида в заданный интервал.
    """
    length = len(sequence)
    if isinstance(length_bounds, (int, float)):
        return length <= length_bounds
    else:
        return length_bounds[0] <= length <= length_bounds[1]

def is_quality_valid(quality, quality_threshold):
    """
    Проверяет среднее качество рида (шкала phred33).
    """
    if not quality:
        return False
    qual_scores = [ord(char) - 33 for char in quality]
    avg_quality = sum(qual_scores) / len(qual_scores)
    return avg_quality >= quality_threshold
