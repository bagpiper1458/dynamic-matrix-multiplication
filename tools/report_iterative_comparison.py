"""Three-way brute force / recursive / iterative comparison from saved samples."""
import argparse
import base64
import csv
import html
import json
import os
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis/iterative"
KINDS = ("toeplitz", "circulant", "unit-circulant")
PHASES = ("baseline", "recursive_prepare", "recursive_apply", "recursive_end", "iterative_layout",
          "iterative_prepare", "iterative_apply", "iterative_hot_end", "iterative_cold_end")


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def html_table(headers, rows):
    return "<table><thead><tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + \
           "</tr></thead><tbody>" + "".join("<tr>" + "".join(f"<td>{html.escape(str(x))}</td>" for x in row) +
                                         "</tr>" for row in rows) + "</tbody></table>"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plots", action="store_true", help="requires matplotlib; also writes standalone HTML")
    args = parser.parse_args()
    samples = defaultdict(list)
    for kind in KINDS:
        with (OUT / f"{kind}-raw.csv").open(encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f):
                samples[(kind, int(row["N"]), row["phase"])].append(float(row["ms_per_call"]))
    medians = {key: statistics.median(values) for key, values in samples.items()}
    with (OUT / "summary.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["kind", "N", "phase", "median_ms", "p25_ms", "p75_ms", "sample_count"])
        for key, values in sorted(samples.items()):
            quartiles = statistics.quantiles(values, method="inclusive") if len(values) > 1 else [values[0]]*3
            writer.writerow([*key, medians[key], quartiles[0], quartiles[2], len(values)])
    with (OUT / "memory.csv").open(encoding="utf-8-sig", newline="") as f:
        memory = {(row["kind"], int(row["N"])): {k: int(v) for k, v in row.items() if k != "kind"}
                  for row in csv.DictReader(f)}
    sizes = sorted({key[1] for key in medians})
    largest = max(sizes)
    three_headers = [f"N={largest}", "比較情境", "Brute force", "Recursive", "Iterative"]
    three_rows = []
    memory_three_rows = []
    for kind in KINDS:
        for label, rec_phase, it_phase in (
            ("矩陣已準備：ms／次", "recursive_apply", "iterative_apply"),
            ("只共用尺寸排程：ms／次", "recursive_end", "iterative_hot_end"),
            ("首次完整計算：ms／次", "recursive_end", "iterative_cold_end"),
        ):
            three_rows.append([kind, label, *(f"{medians[kind, largest, phase]:.6f}"
                                             for phase in ("baseline", rec_phase, it_phase))])
        m = memory[kind, largest]
        for label, rec_key, it_key in (
            ("首次完整計算峰值：bytes", "recursive_cold_peak_bytes", "iterative_cold_peak_bytes"),
            ("已準備，保留 layout 的總量：bytes", "recursive_total_apply_bytes", "iterative_total_apply_with_layout_bytes"),
            ("已準備，釋放 layout 的總量：bytes", "recursive_total_apply_bytes", "iterative_total_apply_without_layout_bytes"),
            ("每次額外空間，含輸出：bytes", "recursive_apply_extra_bytes", "iterative_apply_extra_bytes"),
        ):
            memory_three_rows.append([kind, label, *(f"{m[key]:,}" for key in ("brute_force_peak_bytes", rec_key, it_key))])
    timing_headers = [f"N={largest}，毫秒", *KINDS]
    timing_rows = [[phase, *(f"{medians[k, largest, phase]:.6f}" for k in KINDS)] for phase in PHASES]
    memory_keys = [
        ("Brute force：運算峰值（含輸出）", "brute_force_peak_bytes"),
        ("遞迴：首次完整計算峰值", "recursive_cold_peak_bytes"),
        ("迭代：首次完整計算峰值", "iterative_cold_peak_bytes"),
        ("遞迴：保存的 plan", "recursive_plan_bytes"),
        ("迭代：保存的 plan", "iterative_plan_bytes"),
        ("迭代：可共用的 layout", "layout_retained_bytes"),
        ("迭代：建立 layout 峰值", "layout_build_peak_bytes"),
        ("遞迴：每次額外空間（含輸出）", "recursive_apply_extra_bytes"),
        ("迭代：每次額外空間（含輸出）", "iterative_apply_extra_bytes"),
        ("遞迴：plan＋運算空間", "recursive_total_apply_bytes"),
        ("迭代：plan＋運算空間（已釋放 layout）", "iterative_total_apply_without_layout_bytes"),
        ("迭代：layout＋plan＋運算空間", "iterative_total_apply_with_layout_bytes"),
    ]
    memory_headers = [f"N={largest}，bytes", *KINDS]
    memory_rows = [[label, *(f"{memory[k, largest][key]:,}" for k in KINDS)] for label, key in memory_keys]
    conclusions = []
    for kind in KINDS:
        speedup = medians[kind, largest, "recursive_apply"]/medians[kind, largest, "iterative_apply"]
        m = memory[kind, largest]
        ratio = m["recursive_total_apply_bytes"]/m["iterative_total_apply_with_layout_bytes"]
        bf_speedup = medians[kind, largest, "baseline"]/medians[kind, largest, "iterative_apply"]
        conclusions.append(f"{kind}：矩陣已準備時，iterative 相對 brute force 的加速比為 {bf_speedup:.2f}，相對 recursive 為 {speedup:.2f}；保留 layout 時，recursive／iterative 總 payload 比為 {ratio:.2f}。首次成本另列。")
    metadata = json.loads((OUT / "metadata.json").read_text(encoding="utf-8"))
    method = (
        "Brute force 是直接雙層迴圈逐項計算（原始 CSV 的 baseline），直接讀取結構化矩陣輸入；三版共用相同 scalar kernel，未另加乘法查表或手寫 SIMD。"
        "同一尺寸內使用相同矩陣、向量與 GF(2⁸) scalar kernel（0x11B），先核對正確性再計時。"
        "每項先暖機，隨機排列量測順序；表格與曲線是每輪平均時間的中位數，原始樣本與四分位數另存 CSV。"
        "這是本機固定輸入重複執行的結果，不代表所有硬體或串流工作負載。計時未啟用運算計數器或記憶體 profiler。"
    )
    boundaries = (
        "recursive_prepare 包含建立與釋放遞迴樹；iterative_layout 包含建立與釋放尺寸排程；"
        "iterative_prepare 使用現成 layout，包含建立與釋放矩陣 plan。apply 使用現成 plan，包含工作緩衝區、輸出配置及 checksum。"
        "recursive_end 包含矩陣預計算與一次 apply；iterative_hot_end 共用 layout，但每次重建矩陣 plan；"
        "iterative_cold_end 每次重新建立 layout、plan，再 apply，並計入釋放成本。各階段中位數不可直接相加當作實測端到端時間。"
    )
    memory_note = (
        "記憶體由獨立程式攔截 new/delete，記錄 requested heap payload；已核對 vector capacity 與節點大小公式。"
        "排除輸入、stack、容器物件本體、allocator metadata、作業系統 RSS；不是整個程序的記憶體。"
        "額外空間包括輸出；總量以仍需保留的資料加上該次 apply 峰值計算。layout 可共用，也可在 plan 準備完成後釋放。"
        "layout build 峰值包含暫存 hash table；沒有把預計算成本隱藏在免費快取中。"
    )
    readme = f"""# C++ Brute force／Recursive／Iterative 三版本比較

本機驗證日期：2026-10-09（Asia/Taipei）。目前 iterative 是論文通式與共用計算概念的等價排程重建，並非 Algorithm 1–4 的逐行翻譯；見 [論文與程式逐項對照](../../docs/paper-correspondence.zh-TW.md)。不宣稱已證明時間或空間最優。

{method}

平台：`{metadata['platform']}`；CPU：`{metadata['processor']}`。編譯參數：`{metadata['compile_flags']}`。完整編譯器資訊、執行參數與來源 SHA-256 見 [metadata.json](metadata.json)。

## 最大尺寸的結果

{table(three_headers, three_rows)}

同一個 brute force 數值出現在三種情境，因為它沒有獨立 layout 或矩陣 plan。這不是三次不同的量測。時間沿用既有同批三版本 benchmark 原始樣本；記憶體已補測 brute force 與首次完整計算的峰值。

{chr(10).join('- ' + c for c in conclusions)}

## 各階段明細

{table(timing_headers, timing_rows)}

{boundaries}

## 記憶體

{table(three_headers, memory_three_rows)}

### 保存資料與準備階段明細

{table(memory_headers, memory_rows)}

{memory_note}

Brute force 每次需要 N² 次乘法、N(N−1) 次 XOR；遞迴與迭代都有 3ⁿ−U 次乘法和 3(3ⁿ−2ⁿ) 次向量端 XOR，其中 U=N 只適用 unit-sum circulant，其他情況 U=0。此迭代版每次需要 3ⁿ+N bytes，遞迴版為 3N−2 bytes，brute force 為 N bytes；這些數字都包含輸出。迭代版較小的 plan 和建立排程的峰值需分別評估。

## 圖與原始資料

![各尺寸時間比較](time.png)

![各尺寸記憶體比較](memory.png)

- [可獨立開啟的報告](report.html)
- [時間中位數與四分位數](summary.csv)、[記憶體各階段](memory.csv)
- 原始時間樣本：{', '.join(f'[{k}]({k}-raw.csv)' for k in KINDS)}
- [API、演算法與重現指令](../../docs/iterative.md)

每張子圖都包含三個版本，顏色固定：綠色 brute force、藍色 recursive、橘色 iterative。時間圖三列依序為矩陣已準備、只共用尺寸排程、首次完整計算；記憶體圖依序為首次峰值、已準備且保留排程的總量、每次額外空間。所有圖的 N 與數值軸採對數刻度。
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    if args.plots:
        bundled = ROOT / ".tools/plotting"
        if bundled.is_dir():
            sys.path.insert(0, str(bundled))
        os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".tools/mplconfig"))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                             "axes.spines.right": False, "figure.facecolor": "white"})
        configurations = {
            "time": [
                [("baseline", "Brute force", "#16856c"), ("recursive_apply", "Recursive", "#2463a5"),
                 ("iterative_apply", "Iterative", "#d37020")],
                [("baseline", "Brute force", "#16856c"), ("recursive_end", "Recursive", "#2463a5"),
                 ("iterative_hot_end", "Iterative", "#d37020")],
                [("baseline", "Brute force", "#16856c"), ("recursive_end", "Recursive", "#2463a5"),
                 ("iterative_cold_end", "Iterative", "#d37020")]],
            "memory": [
                [("brute_force_peak_bytes", "Brute force", "#16856c"),
                 ("recursive_cold_peak_bytes", "Recursive", "#2463a5"),
                 ("iterative_cold_peak_bytes", "Iterative", "#d37020")],
                [("brute_force_peak_bytes", "Brute force", "#16856c"),
                 ("recursive_total_apply_bytes", "Recursive", "#2463a5"),
                 ("iterative_total_apply_with_layout_bytes", "Iterative", "#d37020")],
                [("brute_force_peak_bytes", "Brute force", "#16856c"),
                 ("recursive_apply_extra_bytes", "Recursive", "#2463a5"),
                 ("iterative_apply_extra_bytes", "Iterative", "#d37020")]],
        }
        for name, rows in configurations.items():
            fig, axes = plt.subplots(3, 3, figsize=(15, 12), layout="constrained")
            for col, kind in enumerate(KINDS):
                for r, series in enumerate(rows):
                    ax = axes[r, col]
                    for key, label, color in series:
                        values = [medians[kind, n, key] if name == "time" else memory[kind, n][key] for n in sizes]
                        ax.plot(sizes, values, marker="o", markersize=3, linewidth=1.8, color=color, label=label)
                    ax.set_xscale("log", base=2)
                    ax.set_yscale("log")
                    ax.set_xticks([n for n in sizes if n in (1, 4, 16, 64, 256, 1024)])
                    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
                    ax.grid(True, which="major", alpha=.18)
                    ax.set_xlabel("Matrix dimension N")
                    ax.set_ylabel("ms / call (median)" if name == "time" else "Requested heap bytes")
                    scenario = (("Prepared matrix", "Reuse dimension layout", "Cold: include all setup") if name == "time" else
                                ("Cold: include all setup", "Retained + apply (keep layout)", "Extra application buffers"))[r]
                    ax.set_title(kind + "\n" + scenario)
                    ax.legend(fontsize=8, loc="best")
            fig.suptitle("Brute force / Recursive / Iterative | C++, GF(2^8), 0x11B", fontsize=16)
            fig.savefig(OUT / f"{name}.png", dpi=160)
            fig.savefig(OUT / f"{name}.svg")
            plt.close(fig)
        images = "".join(f'<img alt="{name} comparison" src="data:image/png;base64,{base64.b64encode((OUT / (name + ".png")).read_bytes()).decode()}">' for name in ("time", "memory"))
        page = f"""<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Brute force／Recursive／Iterative 三版本比較</title><style>body{{font:16px/1.7 system-ui,sans-serif;max-width:1260px;margin:40px auto;padding:0 24px;color:#17293b;background:#f6f8fb}}h1,h2{{line-height:1.25}}img{{width:100%;background:white;border-radius:12px;margin:18px 0}}table{{border-collapse:collapse;width:100%;background:white;font-size:14px}}td,th{{padding:9px;text-align:left;border-bottom:1px solid #dde3ec}}th{{background:#eaf0f7}}.scroll{{overflow:auto}}.note{{padding:18px;background:#eaf0f7;border-radius:10px}}</style>
<h1>Brute force／Recursive／Iterative 三版本比較</h1><p>2026-10-09 · N=1–{largest} · GF(2⁸), 0x11B</p>
<p>{method}</p><div class="note">{'<br>'.join(conclusions)}</div>
<h2>三個版本，三種使用情境</h2><div class="scroll">{html_table(three_headers, three_rows)}</div>
<p>Brute force 沒有預計算，因此三種情境使用同一量測值。時間沿用同批原始樣本；記憶體已補測 brute force 與首次完整峰值。</p>
{images}<h2>時間：各階段明細</h2><div class="scroll">{html_table(timing_headers, timing_rows)}</div><p>{boundaries}</p>
<h2>三版本記憶體比較</h2><div class="scroll">{html_table(three_headers, memory_three_rows)}</div><p>{memory_note}</p>
<details><summary>保存資料與準備階段明細</summary><div class="scroll">{html_table(memory_headers, memory_rows)}</div></details>
<h2>程式與論文的關係</h2><p>Brute force 逐項計算 y[i] = XOR a[i−j]·b[j]。Recursive 直接使用論文的 F=A₀(B₀⊕B₁)、G=(A₀⊕A₁)B₁、H=(A₀⊕A₂)B₀，再回傳 [F⊕G, F⊕H]。</p>
<p>目前 Iterative 將相同通式的共用係數、向量和、乘積與輸出組合排成迴圈。它不建立論文完整的 N×N SM_A、SM_B 表，而使用 3ⁿ 個不同乘積與共用表。這是等價排程重建，不是 Algorithm 1–4 的逐行翻譯。</p>
<p>新增測試獨立依 Algorithm 1/2 的集合生成步驟建立公式，對 N=1、2、4、8、16、32 做符號消去核對，並核對程式的係數表。<a href="../../docs/paper-correspondence.zh-TW.md">完整推導、4×4 例子與對照表</a>。</p>
<p>資料結構、hash 去重及覆寫排程是重建時的選擇；尚未證明全域最優。</p></html>"""
        (OUT / "report.html").write_text(page, encoding="utf-8")
    print("\n".join(conclusions))


if __name__ == "__main__":
    main()
