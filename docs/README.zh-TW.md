# Dynamic Multiplication 重建專案

這份專案依你提供的論文重新實作，用於展示研究到程式驗證的完整過程。包含 Python 參考版與 C++17 版本；原始碼已遺失，因此會清楚標註為 reconstruction。

## 實作範圍

- GF(2⁸) 上的 Toeplitz／circulant **矩陣乘向量**，尺寸為 2 的次方。
- 依論文三個遞迴子問題的分解，驗證乘法次數從 4ⁿ 降為 3ⁿ。
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

## GitHub 展示

建議 repository 名稱：`dynamic-matrix-multiplication`。Description 可使用：

> Reconstructed Python/C++ implementation of fast Toeplitz and circulant matrix-vector multiplication over GF(2^8), with correctness tests and reproducible benchmarks.

這個專案尚未替你選擇授權條款，也沒有把論文全文、聯絡資料或履歷加入 repository。GitHub 的公開／私人設定與上傳目的地由你指定。
