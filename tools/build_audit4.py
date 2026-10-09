"""Generate reviewable 4x4 trace, formulas and a standalone step-by-step viewer."""
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from dynamic_matrix.audit4 import make_audit


def table(headers, rows):
    return ['| '+' | '.join(headers)+' |', '| '+' | '.join('---' for _ in headers)+' |',
            *['| '+' | '.join(str(cell) for cell in row)+' |' for row in rows]]


def main():
    output=ROOT/'analysis'/'params4'
    output.mkdir(parents=True,exist_ok=True)
    cases=[make_audit(kind) for kind in ['toeplitz','circulant','unit-circulant']]
    (output/'trace.json').write_text(json.dumps(cases,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    lines=['# 4×4 一般解與 PARAMS 對照', '',
        '這是一份可執行、可稽核的排程設計，用來核對論文一般解與 PARAMS 的重用概念。它不是遺失程式碼的逐行還原，也沒有證明全域最小 XOR 次數、最小空間或最佳時間。', '',
        '## 三個計算階段', '',
        '1. **尺寸／矩陣類型固定**：建立 SM_A、SM_B 的符號對照、共用係數的相依關係、九個乘積的索引，以及區塊組合與最後使用時點。沒有任何矩陣或向量數值可在此先算好。',
        '2. **矩陣固定**：用係數 DAG 計算係數組合；只保留乘向量需要的係數。矩陣改變時必須重算，尺寸 schedule 可以保留。',
        '3. **每個新向量**：重新計算 v01、v23、vAll、v13、v02，計算乘積，再組成四個輸出。這些數值不能跨不同向量沿用。', '',
        '## 一般解如何對應', '',
        '使用零起算索引 i=0..3、mask j=0..3，p=i & j：', '',
        '- 係數集合：A(i,j) = XOR { a[p−s] : s 是 j 的子集合位元遮罩 }。',
        '- 向量集合：B(i,3−j) = XOR { b[(~i & j) OR t] : t 是 (3 XOR j) 的子集合位元遮罩 }。',
        '- circulant 額外套用 a[d] = a[d mod 4]。',
        '- y[i] 是這四個 A·B 的 XOR。B 的互補遮罩反映原文「反向配對」；不是把未反向的 A、B 直接逐項相乘。', '',
        '此位元遮罩寫法是本次對論文集合生成規則的明確化。下面以符號多項式消去驗證，避免只靠幾組數值恰好相等。', '',
        '四個輸出的乘積對照：', '',
        *table(['輸出','mask 0','mask 1','mask 2','mask 3'],[
            ['y0','p0','p1','p3','p4'],['y1','p0','p2','p3','p5'],
            ['y2','p0','p1','p6','p7'],['y3','p0','p2','p6','p8']]), '',
        'p0 被四個輸出使用，p1/p2/p3/p6 各被兩個輸出使用，其餘只被一個輸出使用。9 個唯一乘積取代 16 次逐項重算。', '',
        '## 排程與空間口徑', '',
        '本範例順序為：先算 5 個向量中間和，再算 9 個乘積，接著算 6 個半區塊結果與 4 個輸出。每一步先讀取操作數、計算一個 scalar result，將最後一次使用的變數釋放，然後將 result 寫入可重用的 slot。', '',
        '**本排程峰值為 9 個 field-value slots，另需一個暫存 scalar result。** 這不包括輸入向量、持久係數、相依關係、slot 索引和語言物件／容器開銷，也不是 process memory 量測。四個最終輸出包含在 slots 中。沒有宣稱這是最優排程。', '',
        '執行用的 Python trace 保留完整歷史以便檢查，因此實際佔用遠大於 9 bytes。測試另外以九格陣列重播整個排程，刪除死值後再讀下一步，驗證 slot 重用真的可行。', '',
        '矩陣係數表的「保留」是同一矩陣的多個向量可共用；helper 最後使用後可以釋放，但目前生成對照資料時仍保留以利稽核。', '',
        '## 三種案例計數', '',
        *table(['案例','係數預計算 XOR','每個向量 XOR','乘法','identity copy','持久係數值','slot 峰值'],
            [[c['kind'],c['counts']['coefficient_xors'],c['counts']['vector_xors'],c['counts']['multiplications'],c['counts']['identity_copies'],c['counts']['retained_coefficients'],c['counts']['peak_slots']] for c in cases]), '',
        'unit-circulant 的五次係數 XOR 保留與一般 circulant 相同的保守排程；另外驗證列 XOR=1 需三次 XOR，不包含在表中。尚未對已知常數進一步消去係數儲存或預計算。', '',
        'Toeplitz 的 11 次係數 XOR 與舊遞迴版的 12 次不同，因為 a[−1] XOR a[1] 被辨識成共同運算。circulant 還能藉由模 4 等價與共同係數合併到 5 次。此共用表依符號 operand IDs 判斷，不依賴兩個數值偶然相等。', '']
    for case in cases:
        lines += [f"## {case['kind']}", '', '數值以十六進位表示，GF(2⁸) modulus = 0x11B；所有加法均為 XOR。', '', '```text',
                  *[' '.join(f'{v:02X}' for v in row) for row in case['matrix']], '```',
                  'b = '+' '.join(f'{v:02X}' for v in case['vector']),
                  '\n輸出 = '+' '.join(f'{v:02X}' for v in case['output']), '', '### 矩陣係數預計算', '',
                  *table(['步驟','變數','計算','展開集合','值','用途／最後使用'],[
                      [s['step'],s['id'],s['expression'],s['expanded'],f"{s['value']:02X}",s['last_use']] for s in case['matrix_steps']]), '',
                  '持久係數：'+', '.join(f"{v['id']}={v['value']:02X}" for v in case['persistent_coefficients']), '',
                  '### 九個共用乘積', '',
                  *table(['乘積','係數集合','向量集合','值','使用者'],[
                      [p['id'],p['coefficient_expanded'],p['vector_expanded'],f"{p['value']:02X}",', '.join(p['outputs'])] for p in case['products']]), '',
                  '### 一般解的 16 個位置', '',
                  *table(['輸出','mask','A 集合','反向配對 B 集合','共用乘積'],[
                      [t['output'],t['mask'],t['A'],t['B'],t['product']] for t in case['formula_terms']]), '',
                  '### 每個向量的執行與釋放順序', '',
                  *table(['步驟','變數','計算','值','寫入 slot','讀取後釋放','該變數最後使用','步驟後存活數'],[
                      [s['step'],s['id'],s['expression'],f"{s['value']:02X}",s['slot'],', '.join(s['released']) or '—',s['last_use'],len(s['live_after'])] for s in case['vector_steps']]), '',
                  '### 符號驗證：XOR 消去後留下的項', '',
                  *table(['輸出','存活 monomials','結果'],[[s['output'],' XOR '.join(s['surviving_monomials']),'等於對應矩陣列 × b'] for s in case['certificate']]), '',
                  '符號證明使用一般的九個乘積。unit 案例再以列 XOR=1 的條件將四個 full-row-sum 乘法替換成 copy；這個額外條件在數值執行前驗證。', '']
    lines += ['## 重現與驗證', '', '```sh', 'python tools/build_audit4.py', 'python -m unittest discover -s tests -v', '```', '',
              '測試涵蓋：精確的 16→9 公式映射、symbolic monomial cancellation、獨立多項式乘法的 dense oracle、現有迭代草稿對照，以及 150 組隨機案例的九格 register replay。', '',
              '本階段完成的是可稽核設計與 Python 執行模型，尚未將其宣稱為你的原始程式，也未新增 C++ 迭代效能結論。']
    (output/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    template=(ROOT/'tools'/'audit4_template.html').read_text(encoding='utf-8')
    (output/'index.html').write_text(template.replace('__AUDIT_JSON__',json.dumps(cases,ensure_ascii=False).replace('</','<\\/')),encoding='utf-8')
    print('Generated analysis/params4/{index.html,README.md,trace.json}')
    for c in cases: print(c['kind'],c['counts'],c['output'])


if __name__=='__main__': main()
