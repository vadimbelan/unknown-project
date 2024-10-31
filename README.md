**1. Обнаружение ключевых точек**

SIFT (Scale-Invariant Feature Transform)

❌ **2. Сопоставление ключевых точек**

RANSAC

K-ratio test

❌ **3. Построение карт глубины**

Patch Match

Semi-Global Matching

❌ **4. Определение расположения камеры**

Bundle Adjustment

❌ **5. Реконструкция 3D модели**

SfM (Structure from Motion)

## ⚙️ Настройка окружения (MacOS)

```bash
export CPATH=/opt/homebrew/opt/libomp/include:$CPATH
export LIBRARY_PATH=/opt/homebrew/opt/libomp/lib:$LIBRARY_PATH
```

## 🚀 Сборка проекта
```bash
cmake ..
make
```
