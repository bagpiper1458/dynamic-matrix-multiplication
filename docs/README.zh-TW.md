# Dynamic Multiplication 重建專案

這份專案依作者提供的論文實作，展示從研究推導到程式驗證的完整過程，包含 Python 參考版與 C++17 版本。作者後來找回的舊碼用來產生公式並驗證；PARAMS 共用運算當時只有 pseudo code，因此本專案的共用運算與排程屬於依來源新增的實作。

[PARAMS 完整推導](params-derivation.zh-TW.md) 定義 pos、block、預計算與生命週期，並使用舊碼的 lowbit 生成方式核對。

最新的 [三版本效能與記憶體報告](../analysis/packed/report.html) 比較改善後的 C++ brute force、recursive、iterative，並在同次執行中測量原版供對照。原版公式與結果通過驗證；這次改善係數儲存、預計算與內圈成本，不宣稱已證明全域最佳。[新版 API 與取捨](packed.md) 說明共用係數格式、預配置緩衝區，以及較少 XOR 和較低實測時間的差別。

## 實作範圍

- GF(2⁸) 上的 Toeplitz／circulant **矩陣乘向量**，尺寸為 2 的次方。
- 依論文三個遞迴子問題的分解，驗證乘法次數從 4ⁿ 降為 3ⁿ。
- Python 與 C++ 都提供無遞迴的 PARAMS 迭代版：共用尺寸排程、預計算矩陣係數，再以迴圈乘向量。
- circulant 每列元素 XOR 等於 1 時，驗證 3ⁿ−2ⁿ 次乘法的特例。
- 可重複使用的矩陣預處理 plan，分開記錄預處理與每次乘向量的成本。
- 測試、跨語言案例、benchmark、GitHub Actions 設定，以及推導文件。

## 先執行這三步

在專案資料夾開啟終端機，使用已安裝的 Python 3.10 以上版本：

```sh
python -m dynamic_matrix demo
python -m unittest discover -s tests -v
python -m dynamic_matrix benchmark --kind unit-circulant --output benchmarks/local-unit.json
```

如果 `python` 不是可用指令，可改用電腦上的 Python 執行檔完整路徑。C++ 編譯與執行方式請看 [英文 README](../README.md)。

## 面試前值得理解的內容

1. 為什麼 `A0(x0 XOR x1)` 可以讓兩個輸出區塊共用？
2. 為什麼遞迴式是 `M(n)=3M(n−1)`？
3. 列總和為 1 的條件如何傳到 G/H 子問題，又為何不能傳到 F？
4. 為什麼 3(3ⁿ−2ⁿ) 次加法不包含所有前置計算？
5. 為什麼 Python、C++、不同矩陣大小可能有不同的效能交叉點？

完整回答在 [演算法推導](algorithm.md)。原論文的實驗表標示與計數口徑差異記在 [重建紀錄](reconstruction.md)。

[4×4 逐步核對](../analysis/params4/index.html) 展示通式、共用變數與生命週期；[參考迭代版說明](iterative.md) 提供兩種語言的 API；[較早的三版本實測報告](../analysis/iterative/README.md) 保留作為歷史紀錄，不能直接混用不同工作負載的數字。

## 專案與來源

[GitHub repository](https://github.com/bagpiper1458/dynamic-matrix-multiplication) 收錄程式、測試、推導文件、比較方法與原始測量資料。最新比較位於 `analysis/packed`；`analysis/iterative` 和 `analysis/recursive` 保留先前實作的歷史結果。

專案尚未指定散布授權條款。論文的作者資訊與來源記載於英文 README；repository 不包含論文全文或原實驗環境，不宣稱重現論文中的歷史測量數字。
