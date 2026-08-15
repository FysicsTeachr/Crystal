import numpy as np
import time
import os
import shutil

def main():
    current = 0
    shell_counter = 0
    from_char = " "
    dun = 0

    shell = "../shelll.txt"
    preshell = "../preShl.txt"
    measure = "../measur.txt"
    foundp = "../foundp.txt"
    foundling = "foundling.txt"

    while True:
        try:
            with open(measure, 'r') as file:
                lines = file.readlines()
                if len(lines) == 0:
                    time.sleep(1)
                    continue
                current, shell_counter, from_char, dun = lines[0].split()
                current = int(current)
                shell_counter = int(shell_counter)
        except Exception:
            time.sleep(1)
            continue

        if current == 0:
            continue

        with open(measure, 'w') as file:
            if current == shell_counter:
                file.write(f"0 {shell_counter} s {dun}\n")
            else:
                file.write(f"{current + 1} {shell_counter} s {dun}\n")
        break

    with open("status.txt", 'w') as file:
        file.write("r\n")

    qq2 = np.zeros(14, dtype=np.int8)
    try:
        with open(shell, 'r') as file:
            lines = file.readlines()
            for i in range(current - 1):
                lines.pop(0)
                lines.pop(0)
                lines.pop(0)
            for i in range(3):
                qq2[5 * i:5 * (i + 1)] = list(map(int, lines.pop(0).split()))
    except Exception as e:
        print(f"Error reading {shell}: {e}")

    with open("coord.txt", 'w') as file:
        file.write(f"{current} {shell_counter}\n")
        for i in range(0, 14, 5):
            file.write(" ".join(map(str, qq2[i:i + 5])) + "\n")

    pre, qq2d = 0, 0
    try:
        with open(preshell, 'r') as file:
            lines = file.readlines()
            pre, qq2d = map(int, lines[current - 1].split())
    except Exception as e:
        print(f"Error reading {preshell}: {e}")

    shutil.copy(f"../prefolder/{pre}.pun", "pregeom.pun")

    qq2 = np.zeros(14, dtype=np.int8)
    if qq2d > 0:
        qq2[abs(qq2d) - 1] = 1
    else:
        qq2[abs(qq2d) - 1] = -1

    qt = np.array([30.0] * 14)
    qt[1] = qt[3] = qt[5] = qt[7] = 45.0
    qt *= qq2.astype(float)

    with open("coord1.txt", 'w') as file:
        file.write(f"{current} {shell_counter}\n")
        for i in range(0, 14, 5):
            file.write(" ".join(map(str, qt[i:i + 5])) + "\n")

    with open("cancel.txt", 'w') as file:
        file.write("n\n")

if __name__ == "__main__":
    main()

