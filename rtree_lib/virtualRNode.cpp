#include "virtualRNode.h"
#include "rtree.h"
#include "rentry.h"
#include "rnode.h"
#include "header.h"

using namespace std;

VirtualRNode::VirtualRNode()
{

}

VirtualRNode::~VirtualRNode()
{
	for (int i = 0; i<m_usedspace; i++)
		delete m_entry[i];
	delete[] m_entry;
}

int VirtualRNode::copyData(const VirtualRNode* source)
{
	m_pageid = source->m_pageid;
	m_level = source->m_level;
	m_parent = source->m_parent;
	m_usedspace = source->m_usedspace;
	(source->m_isLeaf) ? (m_isLeaf = true) : (m_isLeaf = false);

	m_entry = new RtreeNodeEntry*[m_usedspace];
	memset(m_entry, 0, sizeof(RtreeNodeEntry*)*m_usedspace);
	for (int i = 0; i<m_usedspace; i++)
		m_entry[i] = source->m_entry[i]->clone();

	return 0;
}

int VirtualRNode::copyData(const RtreeNode& source)
{
	m_pageid = source.m_pageid;
	m_level = source.m_level;
	m_parent = source.m_parent;
	m_usedspace = source.m_usedspace;
	(source.isLeaf()) ? (m_isLeaf = true) : (m_isLeaf = false);

	return 0;
}

int VirtualRNode::copyEntries(const RtreeNode& source, int numbers)
{
	m_entry = new RtreeNodeEntry*[numbers];
	memset(m_entry, 0, sizeof(RtreeNodeEntry*)*numbers);
	for (int i = 0; i<numbers; i++)
		m_entry[i] = source.m_entry[i]->clone();

	return 0;
}

int VirtualRNode::insertEntry(const RtreeNodeEntry* source)
{
	m_pageid = source->m_id;
	m_level = 0;
	m_parent = 0;
	m_usedspace = 0;
	m_isLeaf = true;

	m_entry = new RtreeNodeEntry*[2];
	memset(m_entry, 0, sizeof(RtreeNodeEntry*)* 2);

	m_entry[m_usedspace++] = source->clone();

	return 0;
}

int VirtualRNode::displayMBR()
{
	cout << " stored id=" << m_entry[0]->m_id << endl;
	cout << "node x0:" << m_entry[0]->m_hc.getLower()[0];
	cout << " node y0:" << m_entry[0]->m_hc.getUpper()[0];
	cout << " node x1:" << m_entry[0]->m_hc.getLower()[1];
	cout << " node y1:" << m_entry[0]->m_hc.getUpper()[1] << endl;

	return 0;
}

int VirtualRNode::isLeaf()
{
	return m_isLeaf;
}