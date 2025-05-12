#include <iostream>
#include <vector>
#include <memory>
#include <unordered_map>
#include <functional>

// MIL Instruction Types (operators)
enum class MILOp
{
    // Flow Control
    mMIL, mFUNC, mBLOCK, mJUMP, mJUMPT, mJUMPF, mRETURN, mLABEL, mSWITCH,

    // Arithmetic and Logical Expressions
    mADD, mSUB, mMUL, mDIV, mMOD, mAND, mOR, mXOR, mSHL, mSHR,
    mLT, mLE, mGT, mGE, mEQ, mNE,

    // Memory Operations
    mLD, mST,

    // Function Calls
    mCALL, mCON, mARG,

    // Other
    mNOP
};

// Node Types
enum class NodeType
{
    FLOW,  // Flow node types (e.g., mMIL, mFUNC, mBLOCK)
    EXPR   // Expression node types (e.g., mADD, mMUL)
};

// Node structure
struct MILNode
{
    MILOp op;                              // Operator type
    NodeType type;                         // Node type (Flow or Expression)
    std::vector<std::shared_ptr<MILNode>> children;  // Child nodes
    std::unordered_map<std::string, int> attributes; // Attributes (e.g., value)
    int value;                             // Immediate value or result storage
    bool condition;                        // Condition (for mJUMPT/mJUMPF)

    MILNode(MILOp opType, NodeType nodeType) : op(opType), type(nodeType), value(0), condition(false)
    {
    }
};

// MIL Virtual Machine
class MILVirtualMachine
{
private:
    int programCounter;
    std::shared_ptr<MILNode> root;          // Root of the MIL tree
    std::unordered_map<int, int> memory;    // Memory for the VM
    std::unordered_map<int, std::shared_ptr<MILNode>> labels; // Labels for jumps

public:
    MILVirtualMachine() : root(nullptr), programCounter(0)
    {
    }

    void setRoot(std::shared_ptr<MILNode> node)
    {
        root = node;
    }

    void execute()
    {
        if (root)
           executeNode(root);
    }

private:
    void executeNode(std::shared_ptr<MILNode> node)
    {
        switch (node->op)
        {
            case MILOp::mMIL:
                for (auto &child : node->children) executeNode(child);
                break;

            case MILOp::mFUNC:
                std::cout << "Executing function block\n";
                for (auto &child : node->children) executeNode(child);
                break;

            case MILOp::mBLOCK:
                for (auto &child : node->children) executeNode(child);
                break;

            case MILOp::mADD:
                executeExpression(node, [](int a, int b) { return a + b; });
                break;

            case MILOp::mSUB:
                executeExpression(node, [](int a, int b) { return a - b; });
                break;

            case MILOp::mMUL:
                executeExpression(node, [](int a, int b) { return a * b; });
                break;

            case MILOp::mDIV:
                executeExpression(node, [](int a, int b)
                {
                    if (b == 0) throw std::runtime_error("Division by zero");
                    return a / b;
                });
                break;

            case MILOp::mMOD:
                executeExpression(node, [](int a, int b) { return a % b; });
                break;

            case MILOp::mJUMP:
                programCounter = node->value;
                break;

            case MILOp::mJUMPT:
                if (node->condition) programCounter = node->value;
                break;

            case MILOp::mJUMPF:
                if (!node->condition) programCounter = node->value;
                break;

            case MILOp::mLD:
                node->value = memory[node->attributes["addr"]];
                break;

            case MILOp::mST:
                memory[node->attributes["addr"]] = node->value;
                break;

            case MILOp::mCALL:
                std::cout << "Function call\n";
                // Implementation of function call
                break;

            case MILOp::mRETURN:
                std::cout << "Returning from function\n";
                break;

            case MILOp::mSWITCH:
                executeSwitch(node);
                break;

            case MILOp::mNOP:
                break;

            default:
                throw std::runtime_error("Unknown MIL operation");
        }
        programCounter++;
    }

    void executeExpression(std::shared_ptr<MILNode> node, std::function<int(int, int)> operation)
    {
        if (node->children.size() < 2)
            throw std::runtime_error("Invalid expression node");

        auto left = node->children[0];
        auto right = node->children[1];

        executeNode(left);
        executeNode(right);

        node->value = operation(left->value, right->value);
    }

    void executeSwitch(std::shared_ptr<MILNode> node)
    {
        if (node->children.empty())
            throw std::runtime_error("Empty SWITCH block");

        int conditionValue = node->attributes["condition"];
        for (auto &child : node->children)
        {
            if (child->attributes["case"] == conditionValue)
            {
                executeNode(child);
                break;
            }
        }
    }
};

// Example Program Construction
int main()
{
    // Create a MIL program tree
    auto root = std::make_shared<MILNode>(MILOp::mMIL, NodeType::FLOW);

    // Function node
    auto funcNode = std::make_shared<MILNode>(MILOp::mFUNC, NodeType::FLOW);
    root->children.push_back(funcNode);

    // Block node
    auto blockNode = std::make_shared<MILNode>(MILOp::mBLOCK, NodeType::FLOW);
    funcNode->children.push_back(blockNode);

    // Add node
    auto addNode = std::make_shared<MILNode>(MILOp::mADD, NodeType::EXPR);
    addNode->attributes["left"] = 5;  // Left operand
    addNode->attributes["right"] = 3; // Right operand
    blockNode->children.push_back(addNode);

    // Set up the MIL VM
    MILVirtualMachine vm;
    vm.setRoot(root);

    // Execute
    vm.execute();

    return 0;
}

