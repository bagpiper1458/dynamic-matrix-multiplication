"""Plot independently measured recursive timing and memory; no iterative results."""
import base64
import csv
import hashlib
import html
import json
import math
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis' / 'recursive'
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.tools' / 'mplconfig'))
# Optional project-local plotting dependency; never needed by the core library.
if (ROOT / '.tools' / 'plotting').exists():
    sys.path.insert(0, str(ROOT / '.tools' / 'plotting'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker

font = Path('C:/Windows/Fonts/msjh.ttc')
if font.exists():
    font_manager.fontManager.addfont(str(font))
    plt.rcParams['font.family'] = [font_manager.FontProperties(fname=str(font)).get_name(), 'DejaVu Sans']
plt.rcParams.update({'font.size': 11, 'axes.unicode_minus': False, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.facecolor': '#f5f7fa', 'axes.facecolor': '#ffffff',
                     'savefig.facecolor': '#f5f7fa'})


def csv_rows(name):
    with (OUT / name).open(encoding='utf-8-sig', newline='') as f:
        return [{k: (v if k == 'kind' else float(v)) for k, v in r.items()} for r in csv.DictReader(f)]


def setup(ax, title, ylabel):
    ax.set_xscale('log', base=2)
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda value, position: f'{value:g}'))
    ax.set_xticks([1, 4, 16, 64, 256, 1024], ['1', '4', '16', '64', '256', '1024'])
    ax.set_xlabel('矩陣邊長 N（矩陣大小為 N × N）')
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc='left', fontweight='bold', pad=16)
    ax.grid(which='major', color='#dce2eb', linewidth=.7)


def image_uri(name):
    return 'data:image/png;base64,' + base64.b64encode((OUT/name).read_bytes()).decode()


def main():
    py = json.loads((OUT/'python-time.json').read_text(encoding='utf-8'))
    pym = json.loads((OUT/'python-memory.json').read_text(encoding='utf-8'))
    cpp = csv_rows('cpp-time.csv')
    cppm = csv_rows('cpp-memory.csv')
    times = [('C++17 / Clang 18.1.6 / -O3', cpp), ('Python 3.12.14', py['rows'])]
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
    fig.subplots_adjust(left=.07, right=.98, top=.77, bottom=.17, wspace=.25)
    fig.suptitle('遞迴版｜預處理與每次運算的時間分開看', x=.07, y=.97, ha='left', fontsize=22, fontweight='bold')
    fig.text(.07, .88, '一般 circulant 矩陣 · GF(2⁸) · N = 1–1024 · 7 次取中位數 · 對數座標', color='#465367', fontsize=12)
    for ax, (title, rows) in zip(axes, times):
        setup(ax, title, '時間（毫秒／次；越低越好）')
        for key, label, color, style in [
            ('baseline_ms', '傳統法', '#64748b', '--'),
            ('prepare_ms', '只做預處理（含釋放）', '#d97706', ':'),
            ('apply_ms', '遞迴：重用已準備矩陣', '#2563eb', '-'),
            ('end_to_end_ms', '遞迴：預處理＋乘向量', '#dc2626', '-')]:
            ax.plot([r['N'] for r in rows], [r[key] for r in rows], style, color=color, label=label, marker='o', markersize=4, linewidth=2)
        ax.legend(fontsize=9, loc='upper left', framealpha=.92)
    fig.text(.07, .045, '注意：C++ 在這次 N ≤ 1024 的量測中，計入每次預處理仍比傳統法慢。小尺寸的時間受呼叫與計時開銷影響。', fontsize=11, color='#465367')
    fig.savefig(OUT/'recursive-time.png', dpi=160)
    fig.savefig(OUT/'recursive-time.svg')
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
    fig.subplots_adjust(left=.07, right=.98, top=.77, bottom=.19, wspace=.25)
    fig.suptitle('遞迴版｜儲存計算樹，才是主要記憶體成本', x=.07, y=.97, ha='left', fontsize=22, fontweight='bold')
    fig.text(.07, .88, '持久保存的 plan 與單次運算的額外空間分列；不是整個程式的 RSS，也不包含既有輸入。', color='#465367', fontsize=12)
    for ax, title, rows, keys in [
        (axes[0], 'C++：追蹤 new 配置的有效 bytes', cppm,
         ('plan_retained_heap_bytes','prepare_peak_heap_bytes','apply_peak_extra_heap_bytes','baseline_peak_extra_heap_bytes')),
        (axes[1], 'Python：tracemalloc 追蹤配置', pym['rows'],
         ('plan_retained_bytes','prepare_peak_bytes','apply_peak_extra_bytes','baseline_peak_extra_bytes'))]:
        setup(ax, title, '記憶體（KiB；越低越好）')
        for key, label, color, style in zip(keys,
            ('預處理後保留的 plan', '建立 plan 期間的峰值', '乘向量時額外峰值', '傳統法額外峰值'),
            ('#2563eb', '#d97706', '#059669', '#64748b'), ('-', '--', '-', ':')):
            ax.plot([r['N'] for r in rows], [r[key]/1024 for r in rows], style, color=color, label=label, marker='o', markersize=4, linewidth=2)
        ax.legend(fontsize=9, loc='upper left', framealpha=.92)
    fig.text(.07, .055, 'C++ 呼叫深度 n+1：1024×1024 只有 11 層；樹節點卻有 88,573 個。橘線與藍線接近，是因為峰值主要來自樹。', fontsize=11, color='#465367')
    fig.text(.07, .018, 'C++ 排除 allocator 額外管理資訊與呼叫堆疊；Python 與 C++ 的追蹤口徑不同，不能直接拿兩條曲線當作 process memory 比較。', fontsize=10, color='#465367')
    fig.savefig(OUT/'recursive-memory.png', dpi=160)
    fig.savefig(OUT/'recursive-memory.svg')
    plt.close(fig)

    lines = ['# 遞迴版時間與記憶體分析', '', '量測日期：2026-10-07（Asia/Taipei）。本報告僅分析目前重建的遞迴版，沒有混入迭代版，也不代表原論文原始程式的效能。',
        '', '## 主要結論', '',
        '- 時間的理論改善來自三個半尺寸子問題：傳統法 Θ(N²)，預處理後的遞迴運算 Θ(N^log₂3)，約 Θ(N^1.585)。',
        '- N 加倍，理論上傳統乘法量約 4 倍，遞迴乘法量 3 倍；實際時間還包含函式呼叫、索引、配置、快取等成本。',
        '- 目前版本持久儲存完整三叉樹：節點數 (3^(n+1)−1)/2，因此 plan 記憶體為 Θ(N^1.585)。這是實作選擇，不是所有遞迴演算法都必須使用的空間。',
        '- C++ 已使用共用 scratch：單次乘向量為 N-byte 輸出 + (2N−2)-byte 暫存；額外 heap 合計 3N−2 bytes，另有 Θ(log N) 的呼叫堆疊。',
        '- 因此只消除遞迴呼叫，不會自動移除儲存計算樹的成本；迭代版需要同時改變資料表示與共用中間變數的方式。',
        '', '## 時間圖', '', '![時間](recursive-time.png)', '', '## 記憶體圖', '', '![記憶體](recursive-memory.png)',
        '', '## 精確運算與儲存公式', '', '| 項目 | 目前遞迴版，N=2ⁿ | N 加倍時 |', '|---|---|---|',
        '| 一般乘法次數 | 3ⁿ | 3 倍 |',
        '| 列 XOR 為 1 的 circulant 乘法次數 | 3ⁿ−2ⁿ | 趨近 3 倍 |',
        '| 每個向量的加法次數 | 3(3ⁿ−2ⁿ) | 趨近 3 倍 |',
        '| 係數預處理加法次數 | 3·3ⁿ−4·2ⁿ+1 | 趨近 3 倍 |',
        '| 持久樹節點數 | (3^(n+1)−1)/2 | 趨近 3 倍 |',
        '| C++ 單次額外 heap | 3N−2 bytes | 趨近 2 倍 |',
        '| 最深呼叫層數（含 scalar 葉節點） | n+1 | 多 1 層 |',
        '', '列 XOR 為 1 的快捷路徑只減少 N 次 scalar multiplication；目前樹的結構與節點數沒有改變。本次曲線採用一般 circulant，沒有開啟這個快捷路徑。',
        '', '## 實際時間', '', '| 語言 | N | 傳統 ms | 重用 plan ms | 每次重建 ms | 重用加速比 | 每次重建加速比 |', '|---|---:|---:|---:|---:|---:|---:|']
    for label, rows in [('C++',cpp), ('Python',py['rows'])]:
        for r in rows:
            if r['N'] in (16,64,128,256,512,1024):
                lines.append(f"| {label} | {int(r['N'])} | {r['baseline_ms']:.5f} | {r['apply_ms']:.5f} | {r['end_to_end_ms']:.5f} | {r['apply_speedup']:.2f}× | {r['end_to_end_speedup']:.2f}× |")
    lines += ['', '加速比 = 傳統時間 / 遞迴時間；小於 1 表示遞迴版本較慢。這次 C++ 在 N=1024 的 prepared apply 加速約 '
              f"{cpp[-1]['apply_speedup']:.2f}×，但每次重建只有 {cpp[-1]['end_to_end_speedup']:.2f}×。不要把 prepared-only 結果寫成所有使用情境的加速。",
              '', '## 預處理要重用幾次才划算？', '',
              '簡化估算：T_dynamic(K) ≈ P + K·D，T_baseline(K) ≈ K·B。當 B>D 時，嚴格較快的最小整數 K = floor(P/(B−D))+1。P 是含釋放的獨立預處理時間；這是根據分開量測所作的攤提估算，不是 batch benchmark，也未建模 CPU cache/allocator 隨 K 的變化。', '',
              '| 語言 | N=128 | N=256 | N=512 | N=1024 |', '|---|---:|---:|---:|---:|']
    for label, rows in [('C++',cpp), ('Python',py['rows'])]:
        selected = [r for r in rows if r['N'] in (128,256,512,1024)]
        cells = [str(math.floor(r['prepare_ms']/(r['baseline_ms']-r['apply_ms']))+1) if r['baseline_ms']>r['apply_ms'] else '無法攤平' for r in selected]
        lines.append('| ' + ' | '.join([label, *cells]) + ' |')
    lines += ['', '## 實際記憶體', '', '| N | C++ plan MiB | C++ 乘向量額外 bytes | Python plan MiB | Python 乘向量額外 KiB |', '|---:|---:|---:|---:|---:|']
    for c,p in zip(cppm,pym['rows']):
        if c['N'] in (16,64,128,256,512,1024):
            lines.append(f"| {int(c['N'])} | {c['plan_retained_heap_bytes']/2**20:.3f} | {int(c['apply_peak_extra_heap_bytes'])} | {p['plan_retained_bytes']/2**20:.3f} | {p['apply_peak_extra_bytes']/1024:.2f} |")
    lines += ['', '本機 C++ sizeof(Node)=40 bytes。N=1024 時，88,573 個節點即 3,542,920 bytes，約 3.38 MiB；3,070 bytes 是額外 scratch 與輸出，不能誤當作整個演算法的記憶體使用量。傳統 baseline 直接讀取 circulant 第一欄，不建立 N×N dense matrix，所以 baseline 也只需線性空間。',
              '', '## 測量口徑與限制', '',
              f"- OS：{py['platform']}；CPU 識別：{py['processor']}。",
              f"- Python：{py['python']}；C++：Clang 18.1.6 / Zig 0.13.0，`-std=c++17 -O3 -Wall -Wextra -Wpedantic`。",
              '- 時間：每尺寸同一矩陣／向量，先 warmup、隨機交換測量順序，7 次中位數；Python 每次 1 iteration，C++ 每次 100 iterations。兩個語言使用各自的 PRNG；比較的是相同類型與尺寸，沒有宣稱輸入值逐一相同。',
              '- 正確性：每尺寸先比對傳統結果；記憶體測量也另行比對。時間與記憶體追蹤分開執行，避免 profiler 扭曲計時。',
              '- C++ heap：獨立 executable 攔截 new/delete；計 requested bytes，排除 allocator headers、碎片、既有 input、stack、process RSS。追蹤結果另與 node_count × sizeof(Node)、3N−2 公式核對。',
              '- Python：tracemalloc，量測前清除 GC，分階段重新開始追蹤；retained 在清除暫存 free-list 後讀取。排除既有 inputs 與 apply 前已存在的 plan；包含追蹤到的 Python 物件配置，不能等同 C++ 的有效 payload bytes 或 RSS。',
              '- 單次總成本要同時考慮 retained plan + apply extra peak，以及建立 plan 的峰值；不要把不同時間階段的峰值直接全部相加。',
              '- 樹的建立與釋放、遞迴 stack 的常數成本、Python 容器表示，都是此重建版本的特性。尚未逐行還原原論文 PARAMS 版本，亦未證明空間／時間最優。',
              '', '## 後續迭代版的比較基準', '',
              '後續應固定相同輸入與 scalar kernel，分開量測尺寸相關 schedule 建立、矩陣係數 precompute、每個新向量的中間值、單次和多次使用；同時報告共享 schedule、plan、workspace 的記憶體。不能只把 recursion 改成 for-loop，再將差異全部歸因於迭代。',
              '', '## 檔案与重現', '',
              '- `python-time.json` 保留全部時間樣本；CSV 為中位數摘要。C++ CSV 目前只保留中位數，compiler metadata 另存。',
              '- `python-memory.json`、`cpp-memory.csv` 為記憶體原始測量；PNG 與 SVG 由 `tools/plot_recursive_analysis.py` 產出。',
              '- 核心程式與量測程式 SHA-256：', '']
    paths = ['dynamic_matrix/core.py','cpp/dynamic_matrix.hpp','dynamic_matrix/__main__.py','cpp/benchmark.cpp','tools/recursive_memory.py','tools/recursive_memory.cpp']
    for path in paths:
        lines.append(f"  - `{path}`：`{hashlib.sha256((ROOT/path).read_bytes()).hexdigest()}`")
    (OUT/'README.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    # Standalone shareable viewer: both figures embedded, no network dependency.
    summary = f'''<h1>遞迴版：時間與記憶體分析</h1><p>2026-10-07 · 一般 circulant · GF(2⁸) · N = 1–1024</p>
    <p>目前實作以預處理計算樹換取後續運算效率。N 加倍，理論乘法量從傳統的 4 倍降為 3 倍；
    但保存的 plan 也約增為 3 倍。呼叫堆疊僅多一層。</p>
    <img src="{image_uri('recursive-time.png')}" alt="遞迴時間分析圖">
    <p>C++，N=1024：傳統 {cpp[-1]['baseline_ms']:.3f} ms、重用 plan {cpp[-1]['apply_ms']:.3f} ms、每次重建 {cpp[-1]['end_to_end_ms']:.3f} ms。
    重用加速 {cpp[-1]['apply_speedup']:.2f}×；每次重建加速比 {cpp[-1]['end_to_end_speedup']:.2f}×（小於 1 代表較慢）。</p>
    <img src="{image_uri('recursive-memory.png')}" alt="遞迴記憶體分析圖">
    <p>C++，N=1024：plan 約 3.38 MiB，另加 3,070 bytes 輸出／scratch 與 11 層呼叫堆疊。
    這是追蹤配置的量測，不是 RSS；兩個語言的記憶體追蹤口徑也不同。</p>
    <details><summary>完整分析、公式與量測限制</summary><pre>{html.escape(chr(10).join(lines))}</pre></details>'''
    (OUT/'report.html').write_text('<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>遞迴版效能分析</title><style>body{max-width:1200px;margin:36px auto;padding:0 24px;background:#f5f7fa;color:#172337;font:17px/1.7 system-ui}h1{font-size:32px}img{width:100%;height:auto}pre{white-space:pre-wrap;overflow-wrap:anywhere;font:14px/1.7 ui-monospace,monospace}summary{cursor:pointer;font-weight:bold}</style>'+summary+'</html>', encoding='utf-8')
    print('Wrote two PNG/SVG charts, README.md and standalone report.html')


if __name__ == '__main__':
    main()
