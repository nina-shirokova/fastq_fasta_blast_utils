import os


def check_file_exists(file_path: str):
    """
    Проверяет существование файла.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Файл не найден: {file_path}")
    return True


def read_fasta_lines(input_file_path: str):
    """
    Построчное чтение FASTA файла.
    """
    check_file_exists(input_file_path)
    with open(input_file_path, 'r', encoding='utf-8') as file:
        for line in file:
            yield line.strip()


def write_fasta_record(output_file, header: str, sequence: str):
    """
    Записывает заголовок и последовательность в файл.
    """
    output_file.write(f"{header}\n{sequence}\n")


def read_blast_file(blast_file_path: str):
    """
    Читает текстовый файл BLAST и извлекает первые строки значимых выравниваний
    для каждого запроса.
    """
    check_file_exists(blast_file_path)
    hits = []

    with open(blast_file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if "Sequences producing significant alignments:" in line:
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1
            i += 2
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i < len(lines):
                hit_line = lines[i].strip()
                if hit_line:
                    hits.append(hit_line)
        i += 1

    return hits

