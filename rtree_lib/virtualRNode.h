#ifndef VIRTUALNODE_H
#define VIRTUALNODE_H

class RtreeNode;
class RtreeNodeEntry;

class VirtualRNode
{
public:
	int					m_pageid;
	int					m_level;
	int					m_parent;
	int					m_usedspace;
	RtreeNodeEntry**	m_entry;

	bool				m_isLeaf;

	VirtualRNode();
	~VirtualRNode();

	int copyData(const RtreeNode& source);
	int copyData(const VirtualRNode* source);

	int copyEntries(const RtreeNode& source, int numbers);
	int insertEntry(const RtreeNodeEntry* source);
	int isLeaf();

	int displayMBR();
};

#endif