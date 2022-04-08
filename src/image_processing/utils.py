import matplotlib.pyplot as plt
import matplotlib as mpl

FONT_PATH = './fonts/heiti.ttf'

# show a list of images
def show(imgs, col=5, titles=None, gray=True, chinese_title=False):
    
    chinese_font = mpl.font_manager.FontProperties(fname=FONT_PATH)

    row = max((len(imgs) + col - 1) // col, 2)
    plt.figure(figsize=(col * 2.5, row * 2))
    
    for i in range(len(imgs)):
        ax = plt.subplot(row, col, i + 1)
        if gray:
            plt.imshow(imgs[i], 'gray')
        else:
            plt.imshow(imgs[i])
        if titles is not None:
            if chinese_title:
                plt.title(titles[i], fontproperties=chinese_font, fontsize=16)
            else:
                plt.title(titles[i])
        plt.axis("off")