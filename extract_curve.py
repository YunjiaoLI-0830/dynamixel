import os
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev

# ========= settings =========
IMG_PATH = "design1.png"   
OUT_DIR = "./red_out"
EXPECTED_N = 10            
RED_SCORE_THR = 60        
MIN_AREA = 2               
SPLINE_SMOOTH = 0.2        
KAPPA_SAMPLES = 300
DROP_LOWEST_Y_IF_EXTRA = True  
# ============================

os.makedirs(OUT_DIR, exist_ok=True)

def detect_red_by_dominance(bgr):
   
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    r = rgb[:, :, 0].astype(np.int16)
    g = rgb[:, :, 1].astype(np.int16)
    b = rgb[:, :, 2].astype(np.int16)

    score = r - np.maximum(g, b)   # 红色优势
    mask = (score > RED_SCORE_THR).astype(np.uint8) * 255

    # 轻微膨胀让极小红点更容易成连通域
    mask = cv2.dilate(mask, np.ones((3,3), np.uint8), iterations=1)

    nlab, lab, stats, cents = cv2.connectedComponentsWithStats(mask, connectivity=8)

    pts = []
    areas = []
    for i in range(1, nlab):
        area = stats[i, cv2.CC_STAT_AREA]
        if area < MIN_AREA:
            continue
        cx, cy = cents[i]
        pts.append((cx, cy))
        areas.append(area)

    pts = np.array(pts, dtype=float)
    areas = np.array(areas, dtype=float)
    return pts, areas, mask

def keep_expected_points(pts, areas):
    
    if len(pts) <= EXPECTED_N:
        return pts

    # 先按面积降序取前若干（一般足够）
    idx = np.argsort(areas)[::-1]
    pts = pts[idx]
    if len(pts) > EXPECTED_N and DROP_LOWEST_Y_IF_EXTRA:
        # 丢掉 y 最大的（图像坐标 y 向下增大）直到剩 EXPECTED_N
        while len(pts) > EXPECTED_N:
            drop = np.argmax(pts[:,1])
            pts = np.delete(pts, drop, axis=0)
    else:
        pts = pts[:EXPECTED_N]
    return pts

def order_points_nearest_neighbor(pts):
    if len(pts) == 0:
        return pts

    start = int(np.argmax(pts[:,1]))
    used = np.zeros(len(pts), dtype=bool)
    order = [start]
    used[start] = True

    for _ in range(len(pts)-1):
        last = pts[order[-1]]
        cand = np.where(~used)[0]
        d = np.linalg.norm(pts[cand] - last, axis=1)
        nxt = cand[int(np.argmin(d))]
        order.append(nxt)
        used[nxt] = True

    return pts[np.array(order)]

def curvature_from_points(pts_xy):
    x = pts_xy[:,0]
    z = -pts_xy[:,1]

    ds = np.sqrt(np.diff(x)**2 + np.diff(z)**2)
    s = np.concatenate([[0.0], np.cumsum(ds)])
    if s[-1] < 1e-6:
        raise RuntimeError("Degenerate points.")
    u = s / s[-1]

    tck, _ = splprep([x, z], u=u, s=SPLINE_SMOOTH)

    u_dense = np.linspace(0, 1, KAPPA_SAMPLES)
    xd, zd = splev(u_dense, tck)
    x1, z1 = splev(u_dense, tck, der=1)
    x2, z2 = splev(u_dense, tck, der=2)

    kappa = (x1*z2 - z1*x2) / ((x1**2 + z1**2)**1.5 + 1e-9)  # 1/pixel
    return u_dense, kappa, np.column_stack([xd, zd])

def save_debug_overlay(bgr, pts_ord):
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    plt.figure(figsize=(5,4))
    plt.imshow(rgb)
    for i,(x,y) in enumerate(pts_ord):
        plt.scatter([x],[y], s=50)
        plt.text(x+3, y-3, str(i+1), color="yellow", fontsize=12, weight="bold")
    plt.title("Detected red dots (ordered)")
    plt.axis("off")
    out = os.path.join(OUT_DIR, "debug_overlay.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.show()

def main():
    bgr = cv2.imread(IMG_PATH)
    if bgr is None:
        raise FileNotFoundError(IMG_PATH)

    pts, areas, mask = detect_red_by_dominance(bgr)
    print(f"Detected blobs: {len(pts)}")

    # 保存 mask 方便你看分割结果
    cv2.imwrite(os.path.join(OUT_DIR, "mask_red.png"), mask)

    pts = keep_expected_points(pts, areas)
    print(f"Kept points: {len(pts)} (expected {EXPECTED_N})")

    pts_ord = order_points_nearest_neighbor(pts)
    save_debug_overlay(bgr, pts_ord)

    # 导出给 MATLAB：curve.type='points'，用 (x,z)；这里 z=-y
    df_pts = pd.DataFrame({"x": pts_ord[:,0], "z": -pts_ord[:,1]})
    df_pts.to_csv(os.path.join(OUT_DIR, "red_points_ordered.csv"), index=False)

    # 曲率（debug 用）
    u, kappa, fitted = curvature_from_points(pts_ord)
    pd.DataFrame({"s_norm": u, "kappa_1_per_pixel": kappa}).to_csv(
        os.path.join(OUT_DIR, "kappa_profile.csv"), index=False
    )

    plt.figure()
    plt.plot(u, kappa)
    plt.title("Curvature from red dots (spline fit)")
    plt.xlabel("s/L (normalized)")
    plt.ylabel("kappa (1/pixel)")
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(OUT_DIR, "kappa_plot.png"), dpi=200, bbox_inches="tight")
    plt.show()

    print("Saved:")
    print(" - red_out/red_points_ordered.csv   (MATLAB curve.p)")
    print(" - red_out/mask_red.png             (segmentation debug)")
    print(" - red_out/debug_overlay.png        (points + numbering)")
    print(" - red_out/kappa_profile.csv        (kappa(s) debug)")

if __name__ == "__main__":
    main()
