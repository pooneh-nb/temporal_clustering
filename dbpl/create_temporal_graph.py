import xml.dom.minidom


def main():
    doc = xml.dom.minidom.parse("dblp.xml")
    print(doc.nodeName)


if __name__ == '__main__':
    main()
