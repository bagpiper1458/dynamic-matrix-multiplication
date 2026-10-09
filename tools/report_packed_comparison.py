"""Report the compact three-way comparison; never overwrite historical results."""
import argparse
import base64
import csv
import html
import json
import os
from pathlib import Path
import statistics
import sys
from collections import defaultdict
from report_iterative_comparison import table, html_table

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis/packed"
KINDS = ("toeplitz", "circulant", "unit-circulant")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plots", action="store_true")
    args = parser.parse_args()
    samples = defaultdict(list)
    for kind in KINDS:
        with (OUT / f"{kind}-raw.csv").open(encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f):
                samples[kind, int(row["N"]), row["phase"]].append(float(row["ms_per_call"]))
    medians = {key: statistics.median(value) for key, value in samples.items()}
    with (OUT / "summary.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["kind", "N", "phase", "median_ms", "p25_ms", "p75_ms", "sample_count"])
        for key, value in sorted(samples.items()):
            quartiles = statistics.quantiles(value, method="inclusive") if len(value)>1 else [value[0]]*3
            writer.writerow([*key, medians[key], quartiles[0], quartiles[2], len(value)])
    with (OUT / "memory.csv").open(encoding="utf-8-sig", newline="") as f:
        memory = {(r["kind"], int(r["N"])): {k: int(v) for k,v in r.items() if k != "kind"} for r in csv.DictReader(f)}
    sizes = sorted({key[1] for key in medians})
    n = max(sizes)
    headers = [f"N={n}", "情境", "Brute force", "Recursive", "Iterative"]
    timing_rows, memory_rows, comparison_rows = [], [], []
    for kind in KINDS:
        for label, phases in (
            ("已準備，每次配置緩衝區（ms）", ("baseline", "recursive_apply", "iterative_apply")),
            ("已準備，重用緩衝區（ms）", ("baseline_reuse", "recursive_reuse", "iterative_reuse")),
            ("首次：準備＋運算（ms）", ("baseline", "recursive_end", "iterative_end")),
        ):
            timing_rows.append([kind, label, *(f"{medians[kind,n,p]:.6f}" for p in phases)])
        m = memory[kind,n]
        for label, keys in (
            ("首次 heap 峰值（bytes）", ("brute_force_peak_bytes","recursive_cold_peak_bytes","iterative_cold_peak_bytes")),
            ("已準備，保存資料＋運算空間（bytes）", ("brute_force_peak_bytes","recursive_total_apply_bytes","iterative_total_apply_bytes")),
            ("每次額外緩衝區，含輸出（bytes）", ("brute_force_peak_bytes","recursive_apply_extra_bytes","iterative_apply_extra_bytes")),
            ("重用緩衝區後，每次 new 次數", ("brute_reuse_allocations","recursive_reuse_allocations","iterative_reuse_allocations")),
        ):
            memory_rows.append([kind,label,*(f"{m[k]:,}" for k in keys)])
        for method, old, new in (
            ("Brute force apply", "reference_baseline", "baseline"),
            ("Recursive apply", "reference_recursive_apply", "recursive_apply"),
            ("Iterative apply", "reference_iterative_apply", "iterative_apply"),
            ("Recursive 首次", "reference_recursive_end", "recursive_end"),
            ("Iterative 首次", "reference_iterative_cold", "iterative_end"),
        ):
            comparison_rows.append([kind,method,f"{medians[kind,n,old]:.6f}",f"{medians[kind,n,new]:.6f}",
                                    f"{medians[kind,n,old]/medians[kind,n,new]:.2f}"])
    comparison_headers = [f"N={n}", "方法", "原版（ms）", "新版（ms）", "原版／新版時間"]
    metadata = json.loads((OUT / "metadata.json").read_text(encoding="utf-8"))
    method = (
        "三版共用 GF(2⁸) 的 0x11B 八步 scalar kernel、編譯器及 -O3 設定；沒有另加乘法查表或手寫 SIMD，"
        "編譯器可自行向量化。每個尺寸使用相同矩陣與 16 個輪替向量，先驗證所有輸出與係數順序，再暖機、打散測量順序。"
        "表格為每輪平均時間的中位數，原始樣本與四分位數另存。計時不啟用 counts 或 heap profiler。"
    )
    changes = (
        "Recursive 改存連續葉係數，預計算使用遞迴填入預配置陣列，運算仍遞迴呼叫 F/G/H；2×2 底層合併計算。"
        "Iterative 用逐層迴圈直接準備係數，避免建立 hash layout，再依已推導的 PARAMS 三子規則執行。"
        "一般乘積迴圈移除逐元素 identity 判斷；unit 位置由三進位規則產生，內圈使用預先取得的陣列指標。"
        "Brute force 在 circulant 以 bitmask 取代餘數，Toeplitz 將內圈切成兩段。三者都提供配置版與重用緩衝區版。"
    )
    tradeoff = (
        "此次兩種預計算都使用 3·3ⁿ−4N+1 次係數 XOR，沒有做原 hash layout 的係數子式去重。"
        "這是用較規律的數值計算取代大量結構建置的選擇；原版較少的係數 XOR 並未被當作錯誤。"
        "向量端一般情況仍為 3ⁿ 次乘法、3(3ⁿ−N) 次 XOR，unit-sum circulant 為 3ⁿ−N 次乘法。"
    )
    memory_note = (
        "記憶體為獨立 new/delete profiler 的 requested heap payload，排除輸入、stack、容器本體、allocator metadata 與 RSS。"
        "兩種 plan 都只保存 3ⁿ 個係數 bytes；unit 模式不保存額外葉索引表。"
        "Recursive 的工作緩衝區是 N≤2 時 0，其他尺寸 2N−4 bytes；iterative 是 3ⁿ bytes，另加 N bytes 輸出。"
        "重用模式每次零配置不表示零記憶體，其預配置緩衝區已計入保存資料＋運算空間。遞迴還有未列入 heap 的呼叫堆疊。"
    )
    boundary = (
        "已準備的配置版包含輸出／work 配置、運算、checksum 及釋放；重用版包含運算、檢查與 checksum，但緩衝區在計時前配置。"
        "首次包含該方法的係數準備、一次運算與釋放。新版不需要大型尺寸 layout，所以沒有免費隱藏的 hash 預計算。"
        "原版也在同次 benchmark、相同資料上重新測量；歷史報告仍保留，不能直接混用不同工作負載的舊數字。"
    )
    readme = f"""# 經實作改善後的三版本比較

本機驗證：2026-10-09（Asia/Taipei）。原本的公式與結果通過驗證，但原版的指標樹、hash 排程和內圈索引成本並非唯一實作方式。本版改善這些成本；不宣稱已證明全域最優。

{changes}

{tradeoff}

## 時間

{table(headers,timing_rows)}

{method}

{boundary}

## 記憶體

{table(headers,memory_rows)}

{memory_note}

## 同一次執行中的原版對照

{table(comparison_headers,comparison_rows)}

原版／新版時間比大於 1 才表示新版較快。這些比值只代表本機、這些尺寸與本次工作負載；小尺寸計時容易受到配置及計時開銷影響。

## 圖表與重現

![時間](time.png)

![記憶體](memory.png)

平台：`{metadata['platform']}`；CPU：`{metadata['processor']}`。

- [完整 HTML 報告](report.html)、[各階段中位數與四分位數](summary.csv)、[記憶體測量](memory.csv)
- [執行參數、來源 SHA-256 與編譯器資訊](metadata.json)
- [程式/API 與實作取捨](../../docs/packed.md)
- [歷史三版本報告](../iterative/README.md)、[PARAMS 推導](../../docs/params-derivation.zh-TW.md)
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    if args.plots:
        bundled = ROOT / ".tools/plotting"
        if bundled.is_dir(): sys.path.insert(0,str(bundled))
        os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".tools/mplconfig"))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
        plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False})
        specs = {
            "time": [("Prepared: allocate buffers", ("baseline","recursive_apply","iterative_apply")),
                     ("Prepared: reuse buffers", ("baseline_reuse","recursive_reuse","iterative_reuse")),
                     ("Cold: include all preparation", ("baseline","recursive_end","iterative_end"))],
            "memory": [("Cold: peak requested heap", ("brute_force_peak_bytes","recursive_cold_peak_bytes","iterative_cold_peak_bytes")),
                       ("Plan + application buffers", ("brute_force_peak_bytes","recursive_total_apply_bytes","iterative_total_apply_bytes")),
                       ("Extra buffers including output", ("brute_force_peak_bytes","recursive_apply_extra_bytes","iterative_apply_extra_bytes"))],
        }
        colors = ("#16856c","#2463a5","#d37020")
        labels = ("Brute force","Recursive","Iterative")
        for name, rows in specs.items():
            fig, axes = plt.subplots(3,3,figsize=(15,12),layout="constrained")
            for col, kind in enumerate(KINDS):
                for r,(scenario,keys) in enumerate(rows):
                    ax = axes[r,col]
                    for key,color,label in zip(keys,colors,labels):
                        values = [medians[kind,size,key] if name == "time" else memory[kind,size][key] for size in sizes]
                        ax.plot(sizes,values,marker="o",markersize=3,linewidth=1.8,color=color,label=label)
                    ax.set_xscale("log",base=2); ax.set_yscale("log")
                    ax.set_xticks([size for size in sizes if size in (1,4,16,64,256,1024)])
                    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_: f"{x:g}"))
                    ax.grid(True,which="major",alpha=.18)
                    ax.set_title(kind+"\n"+scenario)
                    ax.set_xlabel("Matrix dimension N")
                    ax.set_ylabel("ms / call (median)" if name == "time" else "Requested heap bytes")
                    ax.legend(fontsize=8)
            fig.suptitle("Compact C++ implementations | Brute force / Recursive / Iterative",fontsize=16)
            fig.savefig(OUT / f"{name}.png",dpi=160)
            fig.savefig(OUT / f"{name}.svg")
            plt.close(fig)
        images = "".join(f'<img alt="{name} comparison" src="data:image/png;base64,{base64.b64encode((OUT / (name+".png")).read_bytes()).decode()}">' for name in ("time","memory"))
        page = f"""<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>改善後的三版本比較</title>
<style>body{{font:16px/1.7 system-ui,sans-serif;max-width:1260px;margin:40px auto;padding:0 24px;color:#17293b;background:#f6f8fb}}h1,h2{{line-height:1.25}}img{{width:100%;background:white;border-radius:12px;margin:18px 0}}table{{border-collapse:collapse;width:100%;background:white;font-size:14px}}td,th{{padding:9px;text-align:left;border-bottom:1px solid #dde3ec}}th{{background:#eaf0f7}}.scroll{{overflow:auto}}.note{{padding:18px;background:#eaf0f7;border-radius:10px}}</style>
<h1>Brute force／Recursive／Iterative：改善後的比較</h1><p>2026-10-09 · C++17 · GF(2⁸), 0x11B · N=1–{n}</p>
<div class="note">原版計算正確，但配置方式仍可改善。新版比較使用連續係數陣列、相同 scalar kernel，並分開計算準備與配置成本；不宣稱全域最優。</div>
<p>{changes}</p><p>{tradeoff}</p><h2>時間</h2><div class="scroll">{html_table(headers,timing_rows)}</div><p>{method}</p><p>{boundary}</p>
{images}<h2>記憶體</h2><div class="scroll">{html_table(headers,memory_rows)}</div><p>{memory_note}</p>
<h2>原版與新版：同次、同資料量測</h2><div class="scroll">{html_table(comparison_headers,comparison_rows)}</div>
<p>結果只適用本次機器與工作負載。<a href="../../docs/packed.md">API、測試與重現方式</a>；<a href="../../docs/params-derivation.zh-TW.md">PARAMS 推導</a>。</p></html>"""
        (OUT / "report.html").write_text(page,encoding="utf-8")
    print(table(headers,[r for r in timing_rows if r[0] == "circulant"]))


if __name__ == "__main__": main()
